import glob
import os
import platform
import subprocess
import sys
from typing import Iterable

from setuptools import Extension, find_packages, setup

try:
    # at the second run for install we already have numpy installed as denpendency
    from numpy import get_include

    numpy_include = get_include()
    del get_include
except:
    # the first run for egg_info does not use numpy
    numpy_include = ""
    print("WARNING no NumPy found, not for build")

try:
    # same as above
    from Cython.Build import cythonize
except ImportError:

    def cythonize(module_list, **kwargs):
        return list(module_list) if isinstance(module_list, Iterable) else [module_list]

    print("WARNING no Cython found, not for build")

PROFILE = 0  # disabled in release, only use in dev
NATIVE_OPTIM = True  # enabled for better speed
FORCE_CYTHONIZE = True  # force rebuild in release, no need to in dev


def get_version() -> str:
    # ref https://packaging.python.org/guides/single-sourcing-package-version/
    # solution 3
    version = {}
    with open("herdingspikes/version.py", "r") as f:
        exec(f.read(), version)

    try:
        commit = (
            subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=os.path.dirname(__file__)
            )
            .decode("utf8")
            .strip()
        )
    except:
        commit = ""

    if any(cmd in sys.argv for cmd in ("sdist", "bdist", "bdist_wheel")):
        # in dist, include commit hash as file but not in version
        if commit:
            with open("herdingspikes/.commit_version", "w") as f:
                f.write(commit)
        return version["__version__"]
    else:
        # in install, include commit hash in version if possible
        commit = "+git." + commit[:8] if commit else ""
        return version["__version__"] + commit


with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

# Compile C++ code
extra_compile_args = [] * NATIVE_OPTIM
link_extra_args = []

if platform.system() == "Darwin":  # macOS
    extra_compile_args += [
        "-mmacosx-version-min=10.14",
        "-F.",
        "-I/opt/homebrew/opt/libomp/include",
        "-Xclang",
        "-fopenmp",
        "-std=c++17",
        "-O3",
        "-march=native",
        "-mtune=native",
    ]
    link_extra_args += [
        "-stdlib=libc++",
        "-mmacosx-version-min=10.14",
        "-lomp",
        "-L/opt/homebrew/opt/libomp/lib",
        "-Xclang",
        "-fopenmp",
    ]
elif platform.system() == "Windows":  # Windows
    extra_compile_args += ["/std:c++17"]
else:  # Linux
    extra_compile_args += [
        "-fopenmp",
        "-std=c++17",
        "-O3",
        "-march=native",
        "-mtune=native",
    ]
    link_extra_args += ["-fopenmp"]

# lightning detection code
ext_src = ["detect.pyx"]
sources = glob.glob("herdingspikes/detection_lightning/**/[A-Z]*.cpp", recursive=True)
sources += [os.path.join("herdingspikes/detection_lightning", fn) for fn in ext_src]
detect_lightning = Extension(
    name="herdingspikes.detection_lightning.detect",
    sources=sources,
    include_dirs=[numpy_include],
    define_macros=[("CYTHON_TRACE_NOGIL", "1" if PROFILE >= 2 else "0")],
    extra_compile_args=extra_compile_args,
    extra_link_args=link_extra_args,
    language="c++",
)

# original detection code, ship for comatibility
ext_src = ["detect.pyx"]
sources = glob.glob(
    "herdingspikes/detection_localisation/**/[A-Z]*.cpp", recursive=True
)
sources += [os.path.join("herdingspikes/detection_localisation", fn) for fn in ext_src]
detect_ext = Extension(
    name="herdingspikes.detection_localisation.detect",
    sources=sources,
    include_dirs=[numpy_include],
    define_macros=[("CYTHON_TRACE_NOGIL", "1" if PROFILE >= 2 else "0")],
    extra_compile_args=extra_compile_args,
    extra_link_args=link_extra_args,
    language="c++",
)

setup(
    name="herdingspikes",
    version=get_version(),
    description="Efficient spike detection and sorting for dense MEA",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mhhennig/HS2",
    author="Matthias Hennig Lab, University of Edinburgh",
    author_email="m.hennig@ed.ac.uk",
    license="GPLv3",
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: Cython",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3 :: Only",
    ],
    keywords="spikes sorting electrophysiology detection",
    python_requires=">=3.12",
    install_requires=[
        "cython",
        "numpy",
        "scipy",
        "scikit-learn",
        "matplotlib",
        "pandas",
        "tqdm",
        "joblib",
        "spikeinterface[full]",
        "h5py",
    ],
    # extras_require={
    #     'tests': [
    #         'spikeinterface>=0.95',  # TODO: needs to include hs-detection
    #         'requests',
    #         'tqdm',
    #         'gprof2dot',
    #         'flameprof',
    #         'line_profiler',
    #         'py-spy'
    #     ]
    # },
    packages=find_packages(),
    package_data={
        "herdingspikes": [
            ".commit_version",
            "detection_localisation/**",
            "detection_lightning/**",
        ]
    },
    exclude_package_data={
        "herdingspikes": [
            # 'detect/detect.cpp',  # can only be exlcuded from MANIFEST.in
            "detection_localisation/*.so",
            "detection_lightning/*.so",
            "**/__pycache__/*",
        ]
    },
    # zip_safe=False,
    ext_modules=[detect_ext, detect_lightning],
)


try:
    subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=os.path.dirname(__file__))
    # if git success: in git repo, remove file
    os.remove("herdingspikes/.commit_version")
    # if file to remove not exist: still captured by try...except
except:
    # else: keep file, or file not exist
    pass
