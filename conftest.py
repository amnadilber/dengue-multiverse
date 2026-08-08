"""
Test-session setup.

The only job here is import order. The MSVC runtime shim in
``dengue_pk/_msvc_runtime`` has to run before NumPy pins Anaconda's older copy of
``msvcp140.dll``, and under pytest that means before the plugins load — earlier
than any test module or fixture can arrange. A root ``conftest.py`` is the
earliest hook available.

Without this, TensorFlow cannot initialise inside a pytest session even though it
imports perfectly well from a script, and the PINN tests silently skip. A silent
skip is the failure mode worth guarding against: the suite reports success while
testing nothing.
"""

import ctypes
import os
import sys
from pathlib import Path

SRC = Path(__file__).parent / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

# Recorded before the preload, so the test module can report which copy of the
# runtime was already resident when the session started and therefore whether
# the preload could have had any effect.
MSVC_STATE_AT_STARTUP = {"numpy_imported": "numpy" in sys.modules,
                         "resident_path": None}

if os.name == "nt":
    _k32 = ctypes.WinDLL("kernel32", use_last_error=True)
    _k32.GetModuleHandleW.restype = ctypes.c_void_p
    _handle = _k32.GetModuleHandleW("msvcp140.dll")
    if _handle:
        _buf = ctypes.create_unicode_buffer(1024)
        _k32.GetModuleFileNameW(ctypes.c_void_p(_handle), _buf, 1024)
        MSVC_STATE_AT_STARTUP["resident_path"] = _buf.value

    from dengue_pk._msvc_runtime import preload

    preload()

# TensorFlow is imported here, once, at session start.
#
# Importing it from inside a test instead fails on this machine with
# ERROR_DLL_INIT_FAILED, while the identical import succeeds from a plain script
# — pytest's assertion-rewriting import hook interferes with the DLL directory
# handling TensorFlow performs during its own import. Importing before any test
# module is collected sidesteps the hook, and every later `import tensorflow`
# then returns the already-initialised module from sys.modules.
#
# Wrapped because TensorFlow is an optional dependency: without it the PINN
# tests skip, which is correct, and the rest of the suite is unaffected.
try:  # noqa: SIM105
    import tensorflow  # noqa: F401
except Exception:  # pragma: no cover - environment-dependent
    pass
