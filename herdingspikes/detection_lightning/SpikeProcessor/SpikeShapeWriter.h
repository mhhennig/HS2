#ifndef SPIKESHAPEWRITER_H
#define SPIKESHAPEWRITER_H

#include <string>
#include <fstream>

#include "SpikeProcessor.h"
#include "../RollingArray.h"

namespace HSDetection
{
    class SpikeShapeWriter : public SpikeProcessor
    {
    private:
        std::ofstream spikeFile;
        IntVolt *buffer; // created and released here

        const RollingArray *pTrace; // passed in, should not release here

        IntFrame cutoutStart;
        IntFrame cutoutLen; // cutoutStart + 1 + cutoutEnd, 1 is peak frame

    public:
        SpikeShapeWriter(const std::string &filename, const RollingArray *pTrace,
                         IntFrame cutoutStart, IntFrame cutoutEnd);
        ~SpikeShapeWriter();

        using SpikeProcessor::operator(); // allow call on iterator
        void operator()(Spike *pSpike);
    };

} // namespace HSDetection

#endif
