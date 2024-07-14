# distutils: language=c++
# cython: language_level=3

from libc.stdint cimport int16_t, int32_t
from libcpp cimport bool
from libcpp.string cimport string

cdef extern from "Point.h" namespace "HSDetection":
    cdef cppclass Point:
        float x
        float y

cdef extern from "Spike.h" namespace "HSDetection":
    cdef cppclass Spike:
        int32_t frame
        int32_t channel
        int16_t amplitude
        Point position

cdef extern from "Detection.h" namespace "HSDetection":
    cdef cppclass Detection:
        Detection(int32_t numChannels,
                  int32_t chunkSize,
                  int32_t chunkLeftMargin,
                  bool rescale,
                  const float *scale,
                  const float *offset,
                  bool medianReference,
                  bool averageReference,
                  int32_t spikeDur,
                  int32_t ampAvgDur,
                  float threshold,
                  float minAvgAmp,
                  float maxAHPAmp,
                  const float *channelPositions,
                  float neighborRadius,
                  float innerRadius,
                  int32_t temporalJitter,
                  int32_t riseDur,
                  bool decayFiltering,
                  float decayRatio,
                  bool localize,
                  bool saveShape,
                  string filename,
                  int32_t cutoutStart,
                  int32_t cutoutEnd) except +
        void step(float *traceBuffer, int32_t chunkStart, int32_t chunkLen) except +
        int32_t finish() except +
        const Spike *getResult() except +
