#ifndef QUEUEPROCESSOR_H
#define QUEUEPROCESSOR_H

#include "../SpikeQueue.h"

namespace HSDetection
{
    class QueueProcessor
    {
    public:
        QueueProcessor() {}
        virtual ~QueueProcessor() {}

        // copy constructor deleted to protect possible internals
        QueueProcessor(const QueueProcessor &) = delete;
        // copy assignment deleted to protect possible internals
        QueueProcessor &operator=(const QueueProcessor &) = delete;

        virtual void operator()(SpikeQueue *pQueue) = 0;
    };

} // namespace HSDetection

#endif
