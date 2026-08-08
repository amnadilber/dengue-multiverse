"""A multiverse analysis of climate-attribution inference in dengue models.

The package began as a climate-forced transmission study of Pakistani dengue
surveillance, and the name and the worked case both survive from that. What it
became is a robustness study: every usable outbreak in a global surveillance
compilation, fitted under a full factorial of the analysis choices a published
model would not normally state, asking how often the verdict on climate forcing
changes within a single dataset.

Modules:

* ``models``     — the host--vector and SEIR structures, the integrator, R0
* ``climate``    — transmission forcing and the thermal mortality response
* ``inference``  — datasets, likelihoods, multi-start fitting
* ``robustness`` — what counts as an unstable verdict, and how to measure it
* ``locations``  — the climate point for each reporting unit
* ``multipatch`` — the two-patch heterogeneity test
* ``pinn``       — the neural-network comparison, kept as methods validation
"""

from __future__ import annotations

# Must precede every other import: see the module docstring for why the order
# matters. Importing NumPy first would make this a no-op.
from . import _msvc_runtime  # noqa: F401  (imported for its side effect)

import pathlib

import yaml

__version__ = "0.1.0"

ROOT = pathlib.Path(__file__).resolve().parents[2]


def load_config(path: str | pathlib.Path | None = None) -> dict:
    """Load the single configuration file that governs the whole pipeline.

    Scripts call this rather than hard-coding values, so that any number
    appearing in a result can be traced to the configuration that produced it.
    """
    path = pathlib.Path(path) if path else ROOT / "config" / "config.yaml"
    with open(path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def resolve(cfg: dict, key: str) -> pathlib.Path:
    """Resolve a configured relative path against the repository root."""
    p = ROOT / cfg["paths"][key]
    p.mkdir(parents=True, exist_ok=True)
    return p
