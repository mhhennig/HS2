#include <utility>

#include "SpikeFilterer.h"

using namespace std;

namespace HSDetection
{
    SpikeFilterer::SpikeFilterer(const ProbeLayout *pLayout, IntFrame temporalJitter)
        : pLayout(pLayout), temporalJitter(temporalJitter) {}

    SpikeFilterer::~SpikeFilterer() {}

    void SpikeFilterer::operator()(SpikeQueue *pQueue)
    {
        Spike maxSpike = std::move(*pQueue->begin());
        pQueue->erase(pQueue->begin());

        IntFrame frameBound = maxSpike.frame + temporalJitter;
        IntChannel maxChannel = maxSpike.channel;
        IntVolt maxAmp = maxSpike.amplitude;

        pQueue->remove_if([this, frameBound, maxChannel, maxAmp](const Spike &spike)
                          { return spike.frame <= frameBound &&
                                   pLayout->areNeighbors(spike.channel, maxChannel) &&
                                   spike.amplitude <= maxAmp; });
        pQueue->push_front(std::move(maxSpike));
    }

} // namespace HSDetection
