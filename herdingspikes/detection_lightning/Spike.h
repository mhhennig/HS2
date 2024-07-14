#ifndef SPIKE_H
#define SPIKE_H

#include "Point.h"

namespace HSDetection
{
    class Spike
    {
    public:
        IntFrame frame;
        IntChannel channel;
        IntVolt amplitude;
        Point position; // default to (0,0) if not localized

        Spike(IntFrame frame, IntChannel channel, IntVolt amplitude)
            : frame(frame), channel(channel), amplitude(amplitude), position() {}
        ~Spike() {}
    };

} // namespace HSDetection

#endif
