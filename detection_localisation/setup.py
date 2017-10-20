from Cython.Build import cythonize
from setuptools import setup, Extension
import numpy

long_descr = """`See
<http://github.com/>`_.
"""

setup(
    version='1.0',
    author='Oliver Muthmann, Matthias H Hennig, Albert Puente Encinas',
    license='GPL3',
    description='Efficient spike detection for extracellular recordings.',
    long_description=long_descr,
    url='http://github.com/',
    ext_modules=cythonize(Extension(
           "detect",
           sources=["detect.pyx",
                    "SpkDonline.cpp",
                    "SpikeHandler.cpp",
                    "ProcessSpikes.cpp",
                    "FilterSpikes.cpp",
                    "LocalizeSpikes.cpp"],
           language="c++",
           extra_compile_args=['-std=c++11', '-O3'],
    )),
    include_dirs=[numpy.get_include()], requires=['numpy', 'h5py']
)
