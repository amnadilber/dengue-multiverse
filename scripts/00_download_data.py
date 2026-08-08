"""
Pipeline step 0 — acquire raw data.

Downloads are recorded with their source URL, retrieval date and SHA-256
checksum in ``data/raw/PROVENANCE.json``. Anyone re-running this later can
verify they obtained the same bytes; if an upstream file is silently revised,
the checksum will differ and the discrepancy becomes visible rather than
propagating unnoticed into the results.

Existing files are not re-downloaded. Deleting ``data/raw`` forces a fresh
retrieval.
"""

from __future__ import annotations

import hashlib
import json
import sys
import zipfile
from datetime import date

import requests

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1] / "src"))

from dengue_pk import load_config, resolve  # noqa: E402


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def record(provenance: dict, key: str, url: str, path) -> None:
    provenance[key] = {
        "url": url,
        "retrieved": date.today().isoformat(),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def fetch(url: str, dest, timeout: int = 300) -> bool:
    """Download unless the file is already present. Returns True if downloaded."""
    if dest.exists():
        print(f"  present, skipping: {dest.name}")
        return False
    print(f"  downloading {dest.name} ...")
    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        with open(dest, "wb") as fh:
            for chunk in r.iter_content(1 << 20):
                fh.write(chunk)
    print(f"  saved {dest.stat().st_size / 1e6:.1f} MB")
    return True


def main() -> None:
    cfg = load_config()
    raw = resolve(cfg, "raw")
    prov_path = raw / "PROVENANCE.json"
    provenance = json.loads(prov_path.read_text()) if prov_path.exists() else {}

    # --- OpenDengue case counts -------------------------------------------
    od = cfg["data"]["opendengue"]
    print(f"OpenDengue {od['version']}")
    zip_path = raw / od["filename"]
    fetch(od["url"], zip_path)
    record(provenance, "opendengue_zip", od["url"], zip_path)

    csv_path = raw / od["csv_name"]
    if not csv_path.exists():
        print(f"  extracting {od['csv_name']} ...")
        with zipfile.ZipFile(zip_path) as z:
            z.extract(od["csv_name"], raw)
    print(f"  extracted: {csv_path.stat().st_size / 1e6:.0f} MB")

    # --- NASA POWER daily climate, one series per study window -------------
    np_cfg = cfg["data"]["nasa_power"]
    print("\nNASA POWER daily climate")
    for name, w in cfg["windows"].items():
        pt = w["climate_point"]
        # Pad the window so that lagged rainfall has data to draw on before the
        # first modelled week.
        start = w["start"].replace("-", "")
        end = w["end"].replace("-", "")
        start_padded = str(int(start) - 10000)  # one year earlier
        url = (f"{np_cfg['base_url']}?parameters={','.join(np_cfg['parameters'])}"
               f"&community={np_cfg['community']}"
               f"&longitude={pt['lon']}&latitude={pt['lat']}"
               f"&start={start_padded}&end={end}&format=CSV")
        dest = raw / f"climate_{name}_{pt['name'].lower()}.csv"
        print(f"  {name}: {pt['name']} ({pt['lat']}, {pt['lon']})")
        fetch(url, dest)
        record(provenance, f"climate_{name}", url, dest)

    prov_path.write_text(json.dumps(provenance, indent=2))
    print(f"\nProvenance written: {prov_path.relative_to(prov_path.parents[2])}")


if __name__ == "__main__":
    main()
