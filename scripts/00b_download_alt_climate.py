"""
Pipeline step 0b — climate series for alternative locations.

Each study window is forced by the climate of a single representative city, which
is a real simplification: Pakistan's national series aggregates outbreaks from
Karachi to Peshawar, and no one point speaks for all of them. Step 04 tests how
much that choice matters, and it needs the alternatives downloaded first.

Separate from step 00 because these are only needed for the sensitivity analysis,
and because the download is slow over a constrained connection.
"""

from __future__ import annotations

import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402

import requests  # noqa: E402

ALTERNATIVES = {
    "national_2013": [("Lahore", 31.55, 74.35), ("Karachi", 24.86, 67.01),
                      ("Peshawar", 34.01, 71.58), ("Multan", 30.20, 71.45)],
    "sindh_2021": [("Karachi", 24.86, 67.01), ("Hyderabad", 25.40, 68.37),
                   ("Sukkur", 27.71, 68.83)],
    "kp_2021": [("Peshawar", 34.01, 71.58), ("Abbottabad", 34.15, 73.21),
                ("Bannu", 32.99, 70.60)],
}

cfg = load_config()
raw = resolve(cfg, "raw")
np_cfg = cfg["data"]["nasa_power"]
prov_path = raw / "PROVENANCE.json"
provenance = json.loads(prov_path.read_text()) if prov_path.exists() else {}

for name, window in cfg["windows"].items():
    start = str(int(window["start"].replace("-", "")) - 10000)
    end = window["end"].replace("-", "")
    for city, lat, lon in ALTERNATIVES.get(name, []):
        dest = raw / f"climate_alt_{name}_{city.lower()}.csv"
        if dest.exists():
            print(f"  present: {dest.name}")
            continue
        url = (f"{np_cfg['base_url']}?parameters={','.join(np_cfg['parameters'])}"
               f"&community={np_cfg['community']}&longitude={lon}&latitude={lat}"
               f"&start={start}&end={end}&format=CSV")
        print(f"  downloading {name} / {city} ...")
        try:
            with requests.get(url, stream=True, timeout=300) as r:
                r.raise_for_status()
                with open(dest, "wb") as fh:
                    for chunk in r.iter_content(1 << 20):
                        fh.write(chunk)
        except Exception as exc:
            print(f"    failed: {exc}")
            continue
        provenance[f"climate_alt_{name}_{city.lower()}"] = {
            "url": url, "retrieved": date.today().isoformat(),
            "bytes": dest.stat().st_size}
        print(f"    saved {dest.stat().st_size / 1e3:.0f} kB")

prov_path.write_text(json.dumps(provenance, indent=2))
print(f"\nProvenance updated: {prov_path}")
