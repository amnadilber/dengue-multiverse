"""Pipeline step 27 — a second climate series per location, 1 degree away.

The study's own limitations section names one degree of freedom it measures for
Pakistan and not for the world: a single representative point supplies the
climate for an entire country or province. That is the largest remaining choice
the paper admits to and does not quantify, and a reviewer who has read the
limitations will ask about it first.

It is also the choice an independent 2026 study found matters less than
in-sample correlation suggests: comparing two meteorological stations for dengue
forecasting in Thailand, Khamthong and Phramrung report that stronger marginal
climate--dengue associations did not translate into better out-of-sample
prediction. Whether the same holds for *model selection* — which is a different
question from prediction — is what this makes measurable.

The alternative point is the original shifted one degree of latitude, about 110
km. That distance is deliberately modest: it is well inside the area a national
or provincial case series aggregates over, so a verdict that changes at this
distance is changing for a reason no analyst could defend as substantive. Using a
fixed offset rather than a second named city also keeps the comparison
reproducible and free of a gazetteer that would itself be a set of choices.

Downloads are resumable and skip files already present, so an interrupted run
over a slow connection continues rather than restarting. Network-bound by
design: it is meant to be run alongside a fitting job, not instead of one.
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
from dengue_pk.locations import (CLIMATE_OFFSET_DEG, offset_point,  # noqa: E402
                                 point_for)

import pandas as pd  # noqa: E402
import requests  # noqa: E402

LEAD_IN_DAYS = 400

cfg = load_config()
raw = resolve(cfg, "raw")
tables = resolve(cfg, "tables")
clim_dir = raw / "climate_global"
clim_dir.mkdir(parents=True, exist_ok=True)


def main() -> None:
    inv = pd.read_csv(tables / "12_global_windows.csv",
                      parse_dates=["start", "end"])

    jobs: dict[tuple[str, float, float], dict] = {}
    for _, w in inv.iterrows():
        pt = point_for(w["country"], w["unit"], w["level"])
        if pt is None:
            continue
        name, lat, lon = pt
        key = (name, round(lat, 4), round(lon, 4))
        job = jobs.setdefault(key, {"start": w["start"], "end": w["end"],
                                    "windows": 0})
        job["start"] = min(job["start"], w["start"])
        job["end"] = max(job["end"], w["end"])
        job["windows"] += 1

    print(f"{len(jobs)} locations; downloading an alternative point "
          f"{CLIMATE_OFFSET_DEG} deg of latitude away for each\n")

    np_cfg = cfg["data"]["nasa_power"]
    prov_path = raw / "PROVENANCE.json"
    provenance = json.loads(prov_path.read_text()) if prov_path.exists() else {}

    done = skipped = failed = 0
    t0 = time.time()
    for i, ((name, lat, lon), job) in enumerate(sorted(jobs.items()), 1):
        alt_lat, alt_lon = offset_point(lat, lon)
        slug = (f"{name.lower().replace(' ', '_').replace('/', '_')}"
                f"_alt_{alt_lat}_{alt_lon}")
        dest = clim_dir / f"{slug}.csv"
        if dest.exists() and dest.stat().st_size > 2000:
            skipped += 1
            continue

        start = (job["start"].date()
                 - timedelta(days=LEAD_IN_DAYS)).strftime("%Y%m%d")
        end = job["end"].date().strftime("%Y%m%d")
        url = (f"{np_cfg['base_url']}?parameters={','.join(np_cfg['parameters'])}"
               f"&community={np_cfg['community']}&longitude={alt_lon}"
               f"&latitude={alt_lat}&start={start}&end={end}&format=CSV")

        print(f"[{i:3d}/{len(jobs)}] {name:28s} -> {alt_lat},{alt_lon}  "
              f"({(time.time() - t0) / 60:.1f} min elapsed)", flush=True)
        try:
            with requests.get(url, stream=True, timeout=300) as r:
                r.raise_for_status()
                with open(dest, "wb") as fh:
                    for chunk in r.iter_content(1 << 20):
                        fh.write(chunk)
        except Exception as exc:                    # noqa: BLE001
            print(f"          failed: {str(exc)[:70]}", flush=True)
            failed += 1
            continue
        provenance[f"climate_global/{slug}"] = {
            "url": url, "retrieved": date.today().isoformat(),
            "bytes": dest.stat().st_size,
            "note": f"alternative point, {CLIMATE_OFFSET_DEG} deg latitude from "
                    f"{lat},{lon}"}
        done += 1

    prov_path.write_text(json.dumps(provenance, indent=2))
    print(f"\ndownloaded {done}, already present {skipped}, failed {failed}")
    print(f"total time {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
