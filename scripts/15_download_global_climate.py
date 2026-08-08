"""
Pipeline step 14 — climate series for every window in the global inventory.

One NASA POWER request per distinct location rather than per window: several
windows usually share a reporting unit, and a single request spanning the union
of their dates costs the same as one spanning any of them. Each request begins a
year before the earliest window at that location so the lagged rainfall covariate
is available from the first modelled week without back-filling.

Downloads are resumable. Existing files are left alone, so an interrupted run —
likely, over a slow connection — continues rather than restarting.
"""

from __future__ import annotations

import json
import sys
import time
from datetime import date, timedelta
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402
from dengue_pk.locations import point_for  # noqa: E402

import pandas as pd  # noqa: E402
import requests  # noqa: E402

LEAD_IN_DAYS = 400          # covers the rainfall lag and its smoothing window

cfg = load_config()
raw = resolve(cfg, "raw")
tables = resolve(cfg, "tables")
clim_dir = raw / "climate_global"
clim_dir.mkdir(parents=True, exist_ok=True)

inv = pd.read_csv(tables / "12_global_windows.csv",
                  parse_dates=["start", "end"])

# Group windows by the location they will use, and take the union of their dates.
jobs: dict[tuple[str, float, float], dict] = {}
for _, w in inv.iterrows():
    pt = point_for(w["country"], w["unit"], w["level"])
    if pt is None:
        continue
    name, lat, lon = pt
    key = (name, round(lat, 4), round(lon, 4))
    job = jobs.setdefault(key, {"start": w["start"], "end": w["end"],
                                "windows": 0, "units": set()})
    job["start"] = min(job["start"], w["start"])
    job["end"] = max(job["end"], w["end"])
    job["windows"] += 1
    job["units"].add(f"{w['country']}/{w['unit']}")

print(f"{len(inv)} windows -> {len(jobs)} distinct locations to download\n")

np_cfg = cfg["data"]["nasa_power"]
prov_path = raw / "PROVENANCE.json"
provenance = json.loads(prov_path.read_text()) if prov_path.exists() else {}

done = skipped = failed = 0
t0 = time.time()
for i, ((name, lat, lon), job) in enumerate(sorted(jobs.items()), 1):
    slug = f"{name.lower().replace(' ', '_').replace('/', '_')}_{lat}_{lon}"
    dest = clim_dir / f"{slug}.csv"
    if dest.exists() and dest.stat().st_size > 2000:
        skipped += 1
        continue

    start = (job["start"].date() - timedelta(days=LEAD_IN_DAYS)).strftime("%Y%m%d")
    end = job["end"].date().strftime("%Y%m%d")
    url = (f"{np_cfg['base_url']}?parameters={','.join(np_cfg['parameters'])}"
           f"&community={np_cfg['community']}&longitude={lon}&latitude={lat}"
           f"&start={start}&end={end}&format=CSV")

    elapsed = time.time() - t0
    print(f"[{i:3d}/{len(jobs)}] {name:28s} {job['windows']:2d} window(s)  "
          f"{start}-{end}  ({elapsed / 60:.1f} min elapsed)")
    try:
        with requests.get(url, stream=True, timeout=300) as r:
            r.raise_for_status()
            with open(dest, "wb") as fh:
                for chunk in r.iter_content(1 << 20):
                    fh.write(chunk)
    except Exception as exc:
        print(f"          failed: {str(exc)[:70]}")
        failed += 1
        continue
    provenance[f"climate_global/{slug}"] = {
        "url": url, "retrieved": date.today().isoformat(),
        "bytes": dest.stat().st_size,
        "units": sorted(job["units"])[:6]}
    done += 1

prov_path.write_text(json.dumps(provenance, indent=2))
print(f"\ndownloaded {done}, already present {skipped}, failed {failed}")
print(f"total time {(time.time() - t0) / 60:.1f} min")
print(f"climate files in {clim_dir}")
