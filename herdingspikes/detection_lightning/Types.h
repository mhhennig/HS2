#ifndef TYPES_H
#define TYPES_H

#include <cstdint>

// no standard FP type defined by C/C++, float corresponds to np.single and cython.float

namespace HSDetection
{
    // used in interface
    typedef int32_t IntFrame;   // number of frames
    typedef int32_t IntChannel; // number of channels
    typedef int16_t IntVolt;    // quantized voltage
    typedef float FloatRaw;     // raw trace
    typedef float FloatGeom;    // spatial dimension
    typedef float FloatRatio;   // ratio between values
    typedef int32_t IntResult;  // expected number of spikes

    // used only internally
    typedef int_fast64_t IntCalc; // internal calc, at least 64bit and fast

} // namespace HSDetection

#endif
