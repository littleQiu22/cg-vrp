# Compatible with .pyd that is compiled by mingw
import os
mingw_path = os.getenv("MINGW_DLL_PATH")
if mingw_path and os.path.exists(mingw_path):
    os.add_dll_directory(mingw_path)

from ._rcspp import vrptw_pricing
from . import vrptw

__all__ = [
    "vrptw",
    "vrptw_pricing"
]