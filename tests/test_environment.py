"""
Environment checks.

These do not test the science. They test the assumptions the science relies on
holding true in the environment it is run in, which on this machine turned out to
be a real problem rather than a hypothetical one.
"""

import os
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
HEAVY = ("numpy", "pandas", "scipy", "matplotlib", "tensorflow", "yaml")


@pytest.mark.skipif(os.name != "nt", reason="Windows-specific runtime conflict")
def test_system_msvc_runtime_is_preloaded():
    """Importing the package must pin the system MSVC runtime, not Anaconda's.

    Anaconda ships MSVC runtime 14.29 while the system has 14.51, and TensorFlow
    2.21 fails to initialise against the older one. Windows keeps only the first
    DLL loaded under a given name, so the preload has to happen before NumPy
    arrives — which means before anything else in this package.

    This test would have caught the failure directly: the original shim lived in
    the module that needed TensorFlow, by which point NumPy had already pinned
    the wrong copy.
    """
    import dengue_pk  # noqa: F401  (import order is the thing under test)

    from dengue_pk._msvc_runtime import preload
    assert "msvcp140.dll" in preload(), "system MSVC runtime not loadable"


def test_tensorflow_imports_after_the_package():
    """TensorFlow must import once the package has been imported first.

    On this machine that required tracking down an interaction that had nothing
    to do with the science. TensorFlow imported perfectly from a script and
    failed inside pytest with ERROR_DLL_INIT_FAILED. Moving the runtime preload
    progressively earlier — into the module, into the package `__init__`, into a
    root `conftest.py` — changed nothing, because the culprit ran earlier still:
    the `langsmith` pytest plugin, autoloaded from Anaconda's site-packages,
    which pins the older MSVC runtime before any project code executes.
    Confirmed by bisecting the autoloaded plugin list; `-p no:langsmith_plugin`
    is now set in `pyproject.toml`.

    Skipped rather than failed where TensorFlow is genuinely absent: it is an
    optional dependency, needed only for the PINN comparison. The point of this
    test is that the skip should never be caused by an environment conflict
    masquerading as a missing package.
    """
    import dengue_pk  # noqa: F401

    tf = pytest.importorskip("tensorflow")
    assert tf.__version__
    # A trivial op confirms the native runtime actually initialised, which a
    # successful import alone does not guarantee.
    assert float(tf.reduce_sum(tf.constant([1.0, 2.0]))) == pytest.approx(3.0)


@pytest.mark.parametrize("script", sorted((REPO / "scripts").glob("*.py")),
                         ids=lambda p: p.name)
def test_scripts_import_the_package_before_numpy(script):
    """Every script must import `dengue_pk` before NumPy and friends.

    This is the rule the DLL conflict imposes, and it is invisible: a script with
    the imports in the conventional order runs fine until it reaches TensorFlow,
    then fails with a message about a missing Visual C++ redistributable that has
    nothing to do with the real cause. Reordering imports is easy; diagnosing it
    a second time is not, so the ordering is asserted rather than remembered.
    """
    text = script.read_text(encoding="utf-8")
    package_at = None
    heavy_at = None
    for match in re.finditer(r"^\s*(?:from|import)\s+([A-Za-z_][\w.]*)",
                             text, re.MULTILINE):
        module = match.group(1).split(".")[0]
        if module == "dengue_pk" and package_at is None:
            package_at = match.start()
        if module in HEAVY and heavy_at is None:
            heavy_at = match.start()

    if package_at is None or heavy_at is None:
        pytest.skip(f"{script.name} does not import both")
    assert package_at < heavy_at, (
        f"{script.name} imports a compiled dependency before dengue_pk; move the "
        f"`from dengue_pk import ...` lines above them")
