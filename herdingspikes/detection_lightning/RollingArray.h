#ifndef ROLLINGARRAY_H
#define ROLLINGARRAY_H

#include "Types.h"

namespace HSDetection
{
    class RollingArray
    {
    private:
        IntVolt *arrayBuffer; // created and released here

        IntFrame frameMask; // rolling length will be 2^n and mask is 2^n-1 for bit ops
        IntChannel numChannels;

        static constexpr std::align_val_t memAlign = std::align_val_t(512); // align to 4K/8 to avoid 4K alias

        static constexpr IntFrame getMask(IntFrame x) // get minimum 0...01...1 >= x
        {
            x |= x >> 1;
            x |= x >> 2;
            x |= x >> 4;
            x |= x >> 8;  // if int16_t
            x |= x >> 16; // if int32_t
            // x |= x >> 32; // if int64_t
            return x;
        }

    public:
        // RollingArray(IntFrame rollingLen, IntChannel numChannels)
        //     : frameMask(getMask(rollingLen)), numChannels(numChannels)
        // {
        //     arrayBuffer = new (memAlign) IntVolt[(IntCalc)(frameMask + 1) * numChannels];
        // }
        RollingArray(IntFrame rollingLen, IntChannel numChannels)
            : frameMask(getMask(rollingLen)), numChannels(numChannels)
        {
            // std::cout << "Rollingarray\n";
            arrayBuffer = (IntVolt *) operator new[](sizeof(IntVolt) * (IntCalc)(frameMask + 1) * numChannels,memAlign);
        }

        ~RollingArray() { operator delete[](arrayBuffer, memAlign); }

        // copy constructor deleted to protect buffer
        RollingArray(const RollingArray &) = delete;
        // copy assignment deleted to protect buffer
        RollingArray &operator=(const RollingArray &) = delete;

        const IntVolt *operator[](IntFrame frame) const { return arrayBuffer + ((IntCalc)frame & frameMask) * numChannels; }
        IntVolt *operator[](IntFrame frame) { return arrayBuffer + ((IntCalc)frame & frameMask) * numChannels; }

        const IntVolt &operator()(IntFrame frame, IntChannel channel) const { return (*this)[frame][channel]; }
        IntVolt &operator()(IntFrame frame, IntChannel channel) { return (*this)[frame][channel]; }
    };

} // namespace HSDetection

#endif
