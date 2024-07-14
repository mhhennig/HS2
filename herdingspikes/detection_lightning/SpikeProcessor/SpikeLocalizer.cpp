#include <algorithm>

#include "SpikeLocalizer.h"

using namespace std;

namespace HSDetection
{
    SpikeLocalizer::SpikeLocalizer(const ProbeLayout *pLayout, const RollingArray *pTrace,
                                   const RollingArray *pRef, const RollingArray *pBaseline,
                                   IntFrame temporalJitter, IntFrame riseDur)
        : pLayout(pLayout), pTrace(pTrace), pRef(pRef), pBaseline(pBaseline),
          temporalJitter(temporalJitter), riseDur(riseDur) {}

    SpikeLocalizer::~SpikeLocalizer() {}

    void SpikeLocalizer::operator()(Spike *pSpike)
    {
        const vector<IntChannel> &neighbors = pLayout->getInnerNeighbors(pSpike->channel);
        int numNeighbors = neighbors.size();

        vector<IntCalc> weights(numNeighbors);
        for (int i = 0; i < numNeighbors; i++)
        {
            weights[i] = sumCutout(pSpike->frame, neighbors[i]);
        }

        IntCalc median = getMedian(weights);

        Point sumPoint(0, 0);
        FloatGeom sumWeight = 0;
        for (int i = 0; i < numNeighbors; i++)
        {
            IntCalc weight = weights[i] - median; // correction and threshold on median
            if (weight >= 0)
            {
                sumPoint += (weight + eps) * pLayout->getChannelPosition(neighbors[i]);
                sumWeight += weight + eps;
            }
        }

        pSpike->position = sumPoint / sumWeight;
    }

    IntCalc SpikeLocalizer::sumCutout(IntFrame frame, IntChannel channel) const
    {
        IntVolt baseline = (*pBaseline)[frame - riseDur][channel]; // baseline at the start of event

        IntCalc sum = 0;
        for (IntFrame t = frame - temporalJitter; t <= frame + temporalJitter; t++)
        {
            IntVolt volt = (*pTrace)(t, channel) - baseline - (*pRef)(t, 0);
            if (volt > 0)
            {
                sum += volt;
            }
        }
        return sum;
    }

    IntCalc SpikeLocalizer::getMedian(vector<IntCalc> weights) const // copy param to be modified inside
    {
        vector<IntCalc>::iterator middle = weights.begin() + weights.size() / 2;
        nth_element(weights.begin(), middle, weights.end());
        if (weights.size() % 2 == 0)
        {
            return (*middle + *max_element(weights.begin(), middle)) / 2;
        }
        else
        {
            return *middle;
        }
    }

} // namespace HSDetection
