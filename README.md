# Herding Spikes 2 Lightning 

Fast spike sorting for high density multielectrode arrays

[![PyPI version](https://badge.fury.io/py/herdingspikes.svg)](https://badge.fury.io/py/herdingspikes)
[![Build Status](https://travis-ci.org/mhhennig/HS2.svg?branch=master)](https://travis-ci.org/mhhennig/HS2)

![Spikes](documentation/pictures/spikes.png)

---
**Update July 2024**

This is a new version 0.4.000, which introduces new, fast and performant spike detection code.

The final legacy version is 0.3.104, which introduced compatibility with [SpikeInterface](https://github.com/SpikeInterface/spikeinterface). SpikeInterface wraps many spike sorters, can read almost any file format and contains other useful functionality into a single code base.

This new version still supports the old detection code and can be used to transition to the new code.

---

## Software for high density electrophysiology

This software provides functionality for the detection, localisation and clustering of spike data from dense multielectrode arrays based on the methods described in the following papers:

J.-O. Muthmann, H. Amin, E. Sernagor, A. Maccione, D. Panas, L. Berdondini, U.S. Bhalla, M.H. Hennig MH (2015). [Spike detection for large neural populations using high density multielectrode arrays](http://journal.frontiersin.org/article/10.3389/fninf.2015.00028/abstract). Front. Neuroinform. 9:28. doi: 10.3389/fninf.2015.00028.

G. Hilgen, M. Sorbaro, S. Pirmoradian, J.-O. Muthmann, I. Kepiro, S. Ullo, C. Juarez Ramirez, A. Puente Encinas, A. Maccione, L. Berdondini, V. Murino, D. Sona, F. Cella Zanacchi, E. Sernagor, M.H. Hennig (2016). [Unsupervised spike sorting for large scale, high density multielectrode arrays.](http://www.cell.com/cell-reports/fulltext/S2211-1247(17)30236-X) Cell Reports 18, 2521–2532. bioRxiv: <http://dx.doi.org/10.1101/048645>.

This implementation is highly efficient, spike detection and localisation runs in real time on recordings from 4,096 channels at 20+kHz on a desktop PC. Large recordings with millions of events can be sorted in minutes.

Since we believe publicly funded research code should be free and open, all code is released under GPL-3.0.

### Supported systems  <a name="systems"></a>

- any recording system supported by [SpikeInterface](https://github.com/SpikeInterface/spikeinterface)
- [3Brain](http://3brain.com/) BIOCAM and BIOCAM X (custom implementation only in versions 0.3.XXX)
- note this software will only perform well on high density multielectrode arrays, for example the [Neuropixels probe](https://www.neuropixels.org/) - for other systems, use one of the many other sorters available in [SpikeInterface]([)](https://github.com/SpikeInterface/spikeinterface)

## Quick start <a name="quickstart"></a>

The code has been tested with Python version 3.6 and above. If your system does not have Python pre-installed, the [Anaconda distribution](https://www.anaconda.com/download/) may be used.

### All operative systems - Installation via pip

We suggest you use [Miniconda](https://docs.conda.io/en/latest/miniconda.html), [Anaconda](https://www.anaconda.com/download) or [Mamba](https://github.com/mamba-org/mamba) to set up a working Python system. We also recommend installing the code in a virtual environment.

A pip distribution is available and can be installed as follows:

    pip install herdingspikes

Windows users follow the instructions [here](documentation/windows_install.md).

### From source

The module can automatically be installed, including all dependencies, by cloning this repository:

    git clone https://github.com/mhhennig/HS2.git
    
Then run:
    
    pip install numpy cython
    python setup.py install

## Example code

Example code for the different supported systems is in the folder [notebooks](notebooks). These can be run without installing HS2 system-wide, but requires to run ``python setup.py build_ext --inplace`` in the ``HS2`` directory. Next, run ``jupyter notebook`` and navigate to the directory to try the code. Each notebook will download a short segment of raw data.

Go [here](documentation) for documentation. A worked example for Biocam data is [here](documentation/biocam/BioCam-demo.md).

## Contributors, alphabetical <a name="people"></a>

- [Matthias Hennig](http://homepages.inf.ed.ac.uk/mhennig/index.html): Spike sorting
- [Jano Horvath](https://github.com/JanoHorvath): Parameter optimisation
- [Cole Hurwitz](https://github.com/colehurwitz31): Spike detection, localisation and sorting, C++ code
- [Rickey K. Liang](https://lkct.github.io/): Optimised nad refactored spike detection and localisation (Lightning)
- [Oliver Muthmann](mailto:ollimuh@googlemail.com): Spike detection and localisation
- [Albert Puente Encinas](https://github.com/albertpuente): C++ implementation, optimisation and parallelisation
- [Martino Sorbaro](http://martinosorb.github.io): Spike sorting, class structure and much of the python code
- [Cesar Juarez Ramirez](mailto:cesaripn2@gmail.com): Visualisation
- [Raimon Wintzer](https://github.com/lsIand): GUI and visualisation

## Contact

The herders are based at the School of Informatics, University of Edinburgh. Contact us [here](http://homepages.inf.ed.ac.uk/mhennig/contact/), we are happy to help.   
