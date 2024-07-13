# Installing Herding Spikes on a Mac


If compile the C++ code, [Xcode](https://developer.apple.com/xcode/) is required.

This has been tested and works on MacOS 14.5 - if this fails, try to install the OpenMP library. This can be done with the package manager ``brew``. If you do not have ``brew`` installed, you can do so by following the instructions at
[https://brew.sh/](https://brew.sh/).

Then open a terminal and type

    brew install libomp

and install Herding Spikes with

    pip install numpy cython
    pip install herdingspikes

# Installing Herding Spikes on Windows

## 1. Visual Studio

The C++ code in Herding Spikes requires the Microsoft C++ Build tools. Install them from [https://visualstudio.microsoft.com/visual-cpp-build-tools/](https://visualstudio.microsoft.com/visual-cpp-build-tools/). For a minimal setup, choose ``Desktop development with C++``:

<img src="pictures/vs1.png" width="240" />

and select these packages:

<img src="pictures/vs2.png" width="240" />

## 2. Python and Herding Spikes

Install [Anaconda](https://www.anaconda.com/download/#windows) ands create a Python environment. This can be done with the ``Anaconda Navigator`` per mouse click.

Then opoen a ternminal in the newly created environment and type

    pip install numpy cython
    pip install herdingspikes