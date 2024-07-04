# distutils: language=c++
# cython: language_level=3

from libc.stdint cimport int32_t as _int32_t
from libcpp cimport bool as _bool
from libcpp.vector cimport vector

cimport numpy as np

from .Detection cimport Detection, Spike


cdef inline Detection* newDet(_int32_t numChannels,
                              _int32_t chunkSize,
                              _int32_t chunkLeftMargin,
                              _bool rescale,
                              const float *scale,
                              const float *offset,
                              _bool medianReference,
                              _bool averageReference,
                              _int32_t spikeDur,
                              _int32_t ampAvgDur,
                              float threshold,
                              float minAvgAmp,
                              float maxAHPAmp,
                              const float *channelPositions,
                              float neighborRadius,
                              float innerRadius,
                              _int32_t temporalJitter,
                              _int32_t riseDur,
                              _bool decayFiltering,
                              float decayRatio,
                              _bool localize,
                              _bool saveShape,
                              bytes filename,
                              _int32_t cutoutStart,
                              _int32_t cutoutEnd):
    return new Detection(numChannels,
                         chunkSize,
                         chunkLeftMargin,
                         rescale,
                         scale,
                         offset,
                         medianReference,
                         averageReference,
                         spikeDur,
                         ampAvgDur,
                         threshold,
                         minAvgAmp,
                         maxAHPAmp,
                         channelPositions,
                         neighborRadius,
                         innerRadius,
                         temporalJitter,
                         riseDur,
                         decayFiltering,
                         decayRatio,
                         localize,
                         saveShape,
                         filename,
                         cutoutStart,
                         cutoutEnd)

cdef inline void delDet(Detection* det):
    del det
