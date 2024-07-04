#ifndef SPIKELOCALIZER_H
#define SPIKELOCALIZER_H

#include "SpikeProcessor.h"
#include "../ProbeLayout.h"
#include "../RollingArray.h"

namespace HSDetection
{
    class SpikeLocalizer : public SpikeProcessor
    {
    private:
        const ProbeLayout *pLayout;    // passed in, should not release here
        const RollingArray *pTrace;    // passed in, should not release here
        const RollingArray *pRef;      // passed in, should not release here
        const RollingArray *pBaseline; // passed in, should not release here

        IntFrame temporalJitter;
        IntFrame riseDur;

        static constexpr FloatGeom eps = 1e-12;

        IntCalc sumCutout(IntFrame frame, IntChannel channel) const;
        IntCalc getMedian(std::vector<IntCalc> weights) const; // copy param to be modified inside

    public:
        SpikeLocalizer(const ProbeLayout *pLayout, const RollingArray *pTrace,
                       const RollingArray *pRef, const RollingArray *pBaseline,
                       IntFrame temporalJitter, IntFrame riseDur);
        ~SpikeLocalizer();

        using SpikeProcessor::operator(); // allow call on iterator
        void operator()(Spike *pSpike);
    };

} // namespace HSDetection

#endif
