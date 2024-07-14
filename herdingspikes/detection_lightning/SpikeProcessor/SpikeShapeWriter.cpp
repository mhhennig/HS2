#include "SpikeShapeWriter.h"

using namespace std;

namespace HSDetection
{
    SpikeShapeWriter::SpikeShapeWriter(const string &filename, const RollingArray *pTrace,
                                       IntFrame cutoutStart, IntFrame cutoutEnd)
        : spikeFile(filename, ios::binary | ios::trunc),
          buffer(new IntVolt[cutoutStart + 1 + cutoutEnd]), pTrace(pTrace),
          cutoutStart(cutoutStart), cutoutLen(cutoutStart + 1 + cutoutEnd) {}

    SpikeShapeWriter::~SpikeShapeWriter()
    {
        delete[] buffer;
        spikeFile.close();
    }

    void SpikeShapeWriter::operator()(Spike *pSpike)
    {
        IntFrame cutoutStart = pSpike->frame - this->cutoutStart;
        for (IntFrame t = 0; t < cutoutLen; t++)
        {
            buffer[t] = (*pTrace)(cutoutStart + t, pSpike->channel);
        }

        spikeFile.write((const char *)buffer, cutoutLen * sizeof(IntVolt));
    }

} // namespace HSDetection
