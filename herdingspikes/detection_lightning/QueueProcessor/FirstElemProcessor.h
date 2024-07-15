#ifndef FIRSTELEMPROCESSOR_H
#define FIRSTELEMPROCESSOR_H

#include "QueueProcessor.h"
#include "../SpikeProcessor/SpikeProcessor.h"

namespace HSDetection
{
    class FirstElemProcessor : public QueueProcessor
    {
    private:
        SpikeProcessor *pSpkProc; // passed in, should not release here

    public:
        FirstElemProcessor(SpikeProcessor *pSpkProc) : pSpkProc(pSpkProc) {}
        ~FirstElemProcessor() {}

        void operator()(SpikeQueue *pQueue) { (*pSpkProc)(pQueue->begin()); }
    };

} // namespace HSDetection

#endif
