"""
Environment shim: force the system Microsoft Visual C++ runtime to load first.

This machine has Anaconda's bundled MSVC runtime (14.29) shadowing the system one
(14.51). TensorFlow 2.21 requires the newer version and fails to initialise
against the older, with an unhelpful "DLL initialization routine failed".

Import order is what makes this subtle. Windows will not load two DLLs with the
same base name into one process, so whichever copy of ``msvcp140.dll`` arrives
first wins. NumPy and SciPy pull in Anaconda's copy, so preloading the system
copy *after* importing them has no effect — the old one is already resident, and
TensorFlow binds to it.

The preload must therefore happen before any other extension module is imported,
which is why this runs at the very top of ``dengue_pk/__init__.py`` rather than
inside the module that actually needs TensorFlow.

Harmless everywhere else: on a machine without the conflict, or a non-Windows
one, this does nothing.
"""

from __future__ import annotations

import ctypes
import os

_RUNTIME_DLLS = ("msvcp140.dll", "msvcp140_1.dll", "vcruntime140.dll",
                 "vcruntime140_1.dll")


def preload() -> list[str]:
    """Load the system MSVC runtime. Returns the names successfully loaded."""
    if os.name != "nt":
        return []
    system32 = os.path.join(os.environ.get("SystemRoot", r"C:\Windows"),
                           "System32")
    loaded = []
    for name in _RUNTIME_DLLS:
        path = os.path.join(system32, name)
        if not os.path.exists(path):
            continue
        try:
            ctypes.WinDLL(path)
            loaded.append(name)
        except OSError:
            pass
    return loaded


preload()
