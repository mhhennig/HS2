#ifndef PROBELAYOUT_H
#define PROBELAYOUT_H

#include <vector>

#include "Point.h"

namespace HSDetection
{
    class ProbeLayout
    {
    private:
        std::vector<Point> positions;                  // channel positions, NxXY
        std::vector<std::vector<FloatGeom>> distances; // channel distances, NxN

        std::vector<std::vector<IntChannel>> neighborList;      // adjacency list sorted by channel for neighbors, contains self, Nx?
        std::vector<std::vector<IntChannel>> innerNeighborList; // adjacency list sorted by distance for inner neighbors, contains self, Nx?

        FloatGeom neighborRadius;
        FloatGeom innerRadius;

    public:
        ProbeLayout(IntChannel numChannels, const FloatGeom *channelPositions,
                    FloatGeom neighborRadius, FloatGeom innerRadius);
        ~ProbeLayout();

        // copy constructor deleted to avoid copy of containers
        ProbeLayout(const ProbeLayout &) = delete;
        // copy assignment deleted to avoid copy of containers
        ProbeLayout &operator=(const ProbeLayout &) = delete;

        const Point &getChannelPosition(IntChannel channel) const { return positions[channel]; }

        FloatGeom getChannelDistance(IntChannel channel1, IntChannel channel2) const { return distances[channel1][channel2]; }

        const std::vector<IntChannel> &getNeighbors(IntChannel channel) const { return neighborList[channel]; }
        const std::vector<IntChannel> &getInnerNeighbors(IntChannel channel) const { return innerNeighborList[channel]; }

        bool areNeighbors(IntChannel channel1, IntChannel channel2) const { return getChannelDistance(channel1, channel2) < neighborRadius; }
        bool areInnerNeighbors(IntChannel channel1, IntChannel channel2) const { return getChannelDistance(channel1, channel2) < innerRadius; }
        bool areOuterNeighbors(IntChannel channel1, IntChannel channel2) const { return areNeighbors(channel1, channel2) && !areInnerNeighbors(channel1, channel2); }
    };

} // namespace HSDetection

#endif
