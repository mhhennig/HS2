#ifndef DETECTION_H
#define DETECTION_H

#include <string>

#include "ProbeLayout.h"
#include "TraceWrapper.h"
#include "RollingArray.h"
#include "SpikeQueue.h"

namespace HSDetection
{
    class Detection
    {
    private:
        friend SpikeQueue; // allow access to the whole param set

        // constants
        static constexpr IntVolt initBase = 0;  // initial value of baseline
        static constexpr IntVolt initDev = 400; // initial value of deviation
        static constexpr IntVolt tauBase = 4;   // time constant for baseline update
        static constexpr IntVolt devChange = 1; // changing for deviation update
        static constexpr IntVolt minDev = 200;  // minimum level of deviation

        static constexpr IntCalc thrQuant = 256; // 8bit precision

        static constexpr IntChannel channelAlign = 64 / sizeof(IntVolt); // align IntVolt to 64B (assume wider FloatRaw)

        static constexpr IntChannel alignChannel(IntChannel x) { return (x + (channelAlign - 1)) / channelAlign; }

        // input data
        TraceWrapper traceRaw;      // input trace
        IntChannel numChannels;     // number of probe channels
        IntChannel alignedChannels; // number of slices of aligned channels
        IntFrame chunkSize;         // size of each chunk, only the last chunk can be of a different (smaller) size
        IntFrame chunkLeftMargin;   // margin on the left of each chunk

        // rescaling
        bool rescale;       // whether to scale the input
        FloatRaw *scale;    // scale for rescaling
        FloatRaw *offset;   // offset for rescaling
        RollingArray trace; // rescaled and quantized trace to be used

        // common reference
        bool medianReference;   // whether to use CMR (overrides CAR)
        bool averageReference;  // whether to use CAR
        RollingArray commonRef; // common median/average reference

        // running estimation
        RollingArray runningBaseline;  // running estimation of baseline (33 percentile)
        RollingArray runningDeviation; // running estimation of deviation from baseline

        // detection
        IntFrame *spikeTime; // counter for time since spike peak
        IntVolt *spikeAmp;   // spike peak amplitude
        IntCalc *spikeArea;  // area under spike used for average amplitude, actually integral*fps
        bool *hasAHP;        // flag for AHP existence

        IntFrame spikeDur;  // duration of a spike since peak
        IntFrame ampAvgDur; // duration to average amplitude
        IntCalc threshold;  // threshold to detect spikes, used as multiplier of deviation
        IntCalc minAvgAmp;  // threshold for average amplitude of peak, used as multiplier of deviation
        IntCalc maxAHPAmp;  // threshold for voltage level of AHP, used as multiplier of deviation

        // queue processing
        SpikeQueue *pQueue; // spike queue, must be a pointer to be new-ed later

        ProbeLayout probeLayout; // geometry for probe layout

        std::vector<Spike> result; // detection result, use vector to expand as needed

        IntFrame temporalJitter; // temporal jitter of the time of peak in electrical signal
        IntFrame riseDur;        // duration that a spike rises to peak

        // decay filtering
        bool decayFilter;      // whether to use decay filtering instead of normal one
        FloatRatio decayRatio; // ratio of amplitude to be considered as decayed

        // localization
        bool localize; // whether to turn on localization

        // save shape
        bool saveShape;       // whether to save spike shapes to file
        std::string filename; // filename for saving
        IntFrame cutoutStart; // the start of spike shape cutout
        IntFrame cutoutEnd;   // the end of cutout

    private:
        inline void scaleCast(IntVolt *trace, const FloatRaw *input);
        inline void noscaleCast(IntVolt *trace, const FloatRaw *input);
        inline void commonMedian(IntVolt *ref, const IntVolt *trace,
                                 IntVolt *buffer, IntChannel mid);
        inline void commonAverage(IntVolt *ref, const IntVolt *trace);
        void scaleAndAverage(IntFrame chunkStart, IntFrame chunkLen);
        void castAndCommonref(IntFrame chunkStart, IntFrame chunkLen);
        inline void estimation(IntVolt *baselines, IntVolt *deviations,
                               const IntVolt *trace, const IntVolt *ref,
                               const IntVolt *basePrev, const IntVolt *devPrev,
                               IntChannel alignedStart, IntChannel alignedEnd);
        inline void detection(const IntVolt *trace, const IntVolt *ref,
                              const IntVolt *baselines, const IntVolt *deviations,
                              IntChannel channelStart, IntChannel channelEnd, IntFrame t);
        void estimateAndDetect(IntFrame chunkStart, IntFrame chunkLen);

    public:
        Detection(IntChannel numChannels, IntFrame chunkSize, IntFrame chunkLeftMargin,
                  bool rescale, const FloatRaw *scale, const FloatRaw *offset,
                  bool medianReference, bool averageReference,
                  IntFrame spikeDur, IntFrame ampAvgDur,
                  FloatRatio threshold, FloatRatio minAvgAmp, FloatRatio maxAHPAmp,
                  const FloatGeom *channelPositions, FloatGeom neighborRadius, FloatGeom innerRadius,
                  IntFrame temporalJitter, IntFrame riseDur,
                  bool decayFiltering, FloatRatio decayRatio, bool localize,
                  bool saveShape, std::string filename, IntFrame cutoutStart, IntFrame cutoutEnd);
        ~Detection();

        // copy constructor deleted to protect internals
        Detection(const Detection &) = delete;
        // copy assignment deleted to protect internals
        Detection &operator=(const Detection &) = delete;

        void step(FloatRaw *traceBuffer, IntFrame chunkStart, IntFrame chunkLen);
        IntResult finish();
        const Spike *getResult() const;

    }; // class Detection

} // namespace HSDetection

#endif
