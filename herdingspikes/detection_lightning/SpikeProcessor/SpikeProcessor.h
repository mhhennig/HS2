#ifndef SPIKEPROCESSOR_H
#define SPIKEPROCESSOR_H

#include "../SpikeQueue.h"

namespace HSDetection
{
    class SpikeProcessor
    {
    public:
        SpikeProcessor() {}
        virtual ~SpikeProcessor() {}

        // copy constructor deleted to protect possible internals
        SpikeProcessor(const SpikeProcessor &) = delete;
        // copy assignment deleted to protect possible internals
        SpikeProcessor &operator=(const SpikeProcessor &) = delete;

        void operator()(SpikeQueue::iterator itSpike) { (*this)(&*itSpike); };
        virtual void operator()(Spike *pSpike) = 0;
    };

} // namespace HSDetection

#endif
