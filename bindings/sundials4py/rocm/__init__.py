"""ROCm bindings loaded after PyTorch's HIP runtime."""

import ctypes
import sys
from pathlib import Path

import torch

if torch.__version__ != "2.7.0+rocm6.3" or not torch.version.hip:
    raise ImportError(
        "sundials4py-rocm 7.9.0 requires torch 2.7.0+rocm6.3 "
        "from the PyTorch ROCm 6.3 wheel index"
    )

_torch_hip = (Path(torch.__file__).parent / "lib" / "libamdhip64.so").resolve()
if not _torch_hip.is_file():
    raise ImportError(f"PyTorch's HIP runtime was not found at {_torch_hip}")


def _check_hip_runtime():
    """Reject a process that has loaded a different libamdhip64."""
    try:
        maps = Path("/proc/self/maps").read_text()
    except OSError as exc:
        raise ImportError("Cannot verify the loaded HIP runtime") from exc
    paths = {
        Path(line.split(maxsplit=5)[-1]).resolve()
        for line in maps.splitlines()
        if "libamdhip64.so" in line and len(line.split(maxsplit=5)) == 6
    }
    other = paths - {_torch_hip}
    if other:
        raise ImportError(
            "sundials4py and PyTorch would use different HIP runtimes: "
            + ", ".join(map(str, sorted(other)))
        )


_check_hip_runtime()
_hip_library = ctypes.CDLL(str(_torch_hip), mode=ctypes.RTLD_GLOBAL)

from . import _native  # noqa: E402

_check_hip_runtime()
__version__ = _native.__version__
for _name in ("core", "test", "arkode", "cvodes", "idas", "kinsol"):
    _module = getattr(_native, _name)
    globals()[_name] = _module
    sys.modules[f"{__name__}.{_name}"] = _module
