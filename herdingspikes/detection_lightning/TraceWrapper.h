#ifndef TRACEWRAPPER_H
#define TRACEWRAPPER_H

#include "Types.h"

namespace HSDetection
{
    class TraceWrapper
    {
    private:
        FloatRaw *traceBuffer; // passed in, should not release here

        IntFrame frameOffset; // offset of current chunk
        IntChannel numChannels;
        IntFrame chunkSize; // offset advance by this size with each new chunk

    public:
        TraceWrapper(IntFrame leftMargin, IntChannel numChannels, IntFrame chunkSize)
            : traceBuffer(nullptr), frameOffset(-leftMargin - chunkSize),
              numChannels(numChannels), chunkSize(chunkSize) {}
        ~TraceWrapper() {}

        // should be called to both provide a buffer and advance the offset
        void updateChunk(FloatRaw *traceBuffer) { this->traceBuffer = traceBuffer, frameOffset += chunkSize; }

        const FloatRaw *operator[](IntFrame frame) const { return traceBuffer + (IntCalc)(frame - frameOffset) * numChannels; }
        FloatRaw *operator[](IntFrame frame) { return traceBuffer + (IntCalc)(frame - frameOffset) * numChannels; }

        const FloatRaw &operator()(IntFrame frame, IntChannel channel) const { return (*this)[frame][channel]; }
        FloatRaw &operator()(IntFrame frame, IntChannel channel) { return (*this)[frame][channel]; }
    };

} // namespace HSDetection

#endif
