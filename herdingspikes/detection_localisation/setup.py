from Cython.Build import cythonize
from setuptools import setup, Extension
import numpy
# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
import distutils.sysconfig
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")
# ==================================

sources = ["detect.pyx",
           "SpkDonline.cpp",
           "SpikeHandler.cpp",
           "ProcessSpikes.cpp",
           "FilterSpikes.cpp",
           "LocalizeSpikes.cpp"]

setup(
    ext_modules=cythonize(Extension(name="detect",
                                    sources=sources,
                                    language="c++",
                                    extra_compile_args=['-std=c++11', '-O3'],
                                    include_dirs=[numpy.get_include()]
                                    )),
    include_dirs=[numpy.get_include()],
    requires=['numpy', 'h5py'],
    packages=['detect']
)
