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

long_descr = """`See
<http://github.com/>`_.
"""

FOLDER = "herdingspikes/detection_localisation/"
sources = ["detect.pyx",
           "SpkDonline.cpp",
           "SpikeHandler.cpp",
           "ProcessSpikes.cpp",
           "FilterSpikes.cpp",
           "LocalizeSpikes.cpp"]
for i, s in enumerate(sources):
    sources[i] = FOLDER + s

setup(
    version='1.0',
    author='Oliver Muthmann, Matthias H Hennig, Albert Puente Encinas',
    license='GPL3',
    description='Efficient spike detection for extracellular recordings.',
    long_description=long_descr,
    url='http://github.com/mhhennig/HS2/',
    ext_modules=cythonize(Extension(name="herdingspikes.detection_localisation.detect",
                                    sources=sources,
                                    language="c++",
                                    extra_compile_args=['-std=c++11', '-O3'],
                                    include_dirs=[numpy.get_include()]
                                    )),
    # include_dirs=[numpy.get_include()],
    requires=['numpy', 'h5py'],
    packages=['herdingspikes.detection_localisation.detect']
)
