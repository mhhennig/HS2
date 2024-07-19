import subprocess
import os

__version__ = "0.4.001"
base_version = __version__

__commit__ = ""
try:
    _hsPath = os.path.abspath(os.path.dirname(__file__))
    _commitFile = os.path.join(_hsPath, ".commit_version")
    if os.path.exists(_commitFile):
        with open(_commitFile) as f:
            __commit__ = f.read().strip()
    else:
        __commit__ = (
            subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=_hsPath)
            .decode("utf8")
            .strip()
        )
except:
    pass

if __commit__:
    __version__ = __version__ + "+git." + __commit__[:12]

# del _hsPath, _commitFile, os, subprocess

__all__ = ["__version__", "__commit__", "base_version"]
