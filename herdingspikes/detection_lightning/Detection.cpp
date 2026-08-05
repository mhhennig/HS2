#include <algorithm>
#include <numeric>

#include <omp.h>

#include "Detection.h"

using namespace std;

namespace HSDetection
{
    Detection::Detection(IntChannel numChannels, IntFrame chunkSize, IntFrame chunkLeftMargin,
                         bool rescale, const FloatRaw *scale, const FloatRaw *offset,
                         bool lowpass, bool medianReference, bool averageReference,
                         IntFrame spikeDur, IntFrame ampAvgDur,
                         FloatRatio threshold, FloatRatio minAvgAmp, FloatRatio maxAHPAmp,
                         const FloatGeom *channelPositions, FloatGeom neighborRadius, FloatGeom innerRadius,
                         IntFrame temporalJitter, IntFrame riseDur,
                         bool decayFiltering, FloatRatio decayRatio, bool localize,
                         bool saveShape, string filename, IntFrame cutoutStart, IntFrame cutoutEnd)
        : traceRaw(chunkLeftMargin, numChannels, chunkSize),
          numChannels(numChannels), alignedChannels(alignChannel(numChannels)),
          chunkSize(chunkSize), chunkLeftMargin(chunkLeftMargin), rescale(rescale),
          scale((FloatRaw *) operator new[](sizeof(FloatRaw) * alignedChannels * channelAlign, (std::align_val_t(channelAlign * sizeof(IntVolt))))),
          offset((FloatRaw *) operator new[](sizeof(FloatRaw) * alignedChannels * channelAlign, (std::align_val_t(channelAlign * sizeof(IntVolt))))),
          trace(chunkSize + chunkLeftMargin, alignedChannels * channelAlign),
          lowpass(lowpass),
          medianReference(medianReference), averageReference(averageReference),
          commonRef(chunkSize + chunkLeftMargin, 1),
          runningBaseline(chunkSize + chunkLeftMargin, alignedChannels * channelAlign),
          runningDeviation(chunkSize + chunkLeftMargin, alignedChannels * channelAlign),
          spikeTime(new IntFrame[numChannels]), spikeAmp(new IntVolt[numChannels]),
          spikeArea(new IntCalc[numChannels]), hasAHP(new bool[numChannels]),
          spikeDur(spikeDur), ampAvgDur(ampAvgDur), threshold(threshold * thrQuant),
          minAvgAmp(minAvgAmp * thrQuant), maxAHPAmp(maxAHPAmp * thrQuant),
          probeLayout(numChannels, channelPositions, neighborRadius, innerRadius),
          result(), temporalJitter(temporalJitter), riseDur(riseDur),
          decayFilter(decayFiltering), decayRatio(decayRatio), localize(localize),
          saveShape(saveShape), filename(filename), cutoutStart(cutoutStart), cutoutEnd(cutoutEnd)
    {
        fill_n(this->scale, alignedChannels * channelAlign, (FloatRaw)1);
        fill_n(this->offset, alignedChannels * channelAlign, (FloatRaw)0);
        if (rescale)
        {
            copy_n(scale, numChannels, this->scale);
            copy_n(offset, numChannels, this->offset);
        }

        fill_n(runningBaseline[-1], alignedChannels * channelAlign, initBase);
        fill_n(runningDeviation[-1], alignedChannels * channelAlign, initDev);

        fill_n(spikeTime, numChannels, (IntFrame)-1);

        pQueue = new SpikeQueue(this); // all the params should be ready
    }

    Detection::~Detection()
    {
        delete pQueue;

        delete[] spikeTime;
        delete[] spikeAmp;
        delete[] spikeArea;
        delete[] hasAHP;

        operator delete[](scale, align_val_t(channelAlign * sizeof(IntVolt)));
        operator delete[](offset, align_val_t(channelAlign * sizeof(IntVolt)));
    }

    void Detection::step(FloatRaw *traceBuffer, IntFrame chunkStart, IntFrame chunkLen)
    {
        traceRaw.updateChunk(traceBuffer);

#pragma omp parallel
        {
            castAndCommonref(chunkStart, chunkLen);
#pragma omp barrier
            estimateAndDetect(chunkStart, chunkLen);
        }

        pQueue->process();
    }

    IntResult Detection::finish()
    {
        pQueue->finalize();
        return result.size();
    }

    const Spike *Detection::getResult() const
    {
        return result.data();
    }

    // Turn the raw float chunk into the quantized `trace`, in one pass:
    //
    //     src   = lowpass ? (raw[t] + raw[t-1]) / 2 : raw[t]
    //     trace = (src - spatialDC(src)) * scale + offset
    //
    // Three notes on why this is a single expression rather than a matrix of
    // special cases:
    //
    //  - `rescale` needs no branch. When it is off the constructor leaves
    //    scale=1, offset=0, so the affine form is already the identity.
    //  - `offset` is -temporalDC * scale, i.e. the cast is really
    //    (src - temporalDC - spatialDC) * scale. The two DC removals are the
    //    same operation on the two axes; only the spatial one varies per frame.
    //  - the spatial DC must come off BEFORE scaling. It is uniform across
    //    channels only in the raw domain; after per-channel scaling it is
    //    s_c*dc, and median{s*x} != s_c*median{x} unless all s_c are equal.
    //    That is also why it cannot live in `commonRef`: the correction is
    //    per-channel, `commonRef` holds one value per frame.
    //
    // `commonRef` is therefore never written and keeps its zero-fill, so the
    // "- *ref" terms in estimation(), detection() and SpikeLocalizer are "- 0".
    void Detection::castAndCommonref(IntFrame chunkStart, IntFrame chunkLen)
    {
        int numThreads = omp_get_num_threads();
        int threadNum = omp_get_thread_num();
        // each thread have ceil(chunkLen / numThreads), last has fewer
        IntFrame thChunkLen = (chunkLen + numThreads - 1) / numThreads;
        IntFrame thChunkStart = chunkStart + threadNum * thChunkLen;
        thChunkLen = min(thChunkLen, chunkStart + chunkLen - thChunkStart);

        const IntChannel nAligned = alignedChannels * channelAlign;
        const bool needSpatialDC = medianReference || averageReference;
        const IntChannel mid = numChannels / 2;
        // nth_element reorders its input, so the source values get their own buffer
        FloatRaw *buffer = needSpatialDC ? new FloatRaw[numChannels] : nullptr;

        for (IntFrame t = thChunkStart; t < thChunkStart + thChunkLen; t++)
        {
            const FloatRaw *input = traceRaw[t];
            const FloatRaw *inputPrev = lowpass ? traceRaw[t - 1] : nullptr;
            IntVolt *out = trace[t];

            FloatRaw spatialDC = 0;
            if (needSpatialDC)
            {
                if (lowpass)
                {
                    for (IntChannel i = 0; i < numChannels; i++)
                    {
                        buffer[i] = (input[i] + inputPrev[i]) / 2;
                    }
                }
                else
                {
                    copy_n(input, numChannels, buffer);
                }

                if (medianReference)
                {
                    nth_element(buffer, buffer + mid, buffer + numChannels);
                    spatialDC = buffer[mid];
                }
                else
                {
                    spatialDC = accumulate(buffer, buffer + numChannels, (double)0) / numChannels;
                }
            }

            if (lowpass)
            {
                for (IntChannel i = 0; i < nAligned; i++)
                {
                    out[i] = ((input[i] + inputPrev[i]) / 2 - spatialDC) * scale[i] + offset[i];
                }
            }
            else
            {
                for (IntChannel i = 0; i < nAligned; i++)
                {
                    out[i] = (input[i] - spatialDC) * scale[i] + offset[i];
                }
            }
        }

        delete[] buffer;
    }

    void Detection::estimateAndDetect(IntFrame chunkStart, IntFrame chunkLen)
    {
        int numThreads = omp_get_num_threads();
        int threadNum = omp_get_thread_num();
        // each thread have ceil(alignedChannels / numThreads), last has fewer
        IntChannel thChannels = (alignedChannels + numThreads - 1) / numThreads;
        IntChannel thAlignedStart = threadNum * thChannels;
        IntChannel thAlignedEnd = (threadNum + 1) * thChannels;
        thAlignedEnd = min(thAlignedEnd, alignedChannels);
        IntChannel thActualEnd = min(thAlignedEnd * channelAlign, numChannels);

        for (IntFrame t = chunkStart; t < chunkStart + chunkLen; t++)
        {
            estimation(runningBaseline[t], runningDeviation[t],
                       trace[t], commonRef[t],
                       runningBaseline[t - 1], runningDeviation[t - 1],
                       thAlignedStart, thAlignedEnd);

            detection(trace[t], commonRef[t],
                      runningBaseline[t], runningDeviation[t],
                      thAlignedStart * channelAlign, thActualEnd, t);
        }
    }

    void Detection::estimation(IntVolt *baselines, IntVolt *deviations,
                               const IntVolt *trace, const IntVolt *ref,
                               const IntVolt *basePrev, const IntVolt *devPrev,
                               IntChannel alignedStart, IntChannel alignedEnd)
    {
        for (IntChannel i = alignedStart * channelAlign; i < alignedEnd * channelAlign; i++)
        {
            IntVolt volt = trace[i] - *ref - basePrev[i];

            IntVolt dltBase = 0;
            dltBase = (devPrev[i] < volt) ? devPrev[i] / tauBase : dltBase;
            dltBase = (volt < -devPrev[i]) ? -devPrev[i] / (tauBase * 2) : dltBase;
            baselines[i] = basePrev[i] + dltBase;

            IntVolt dltDev = 0;
            dltDev = (devPrev[i] < volt && volt < 5 * devPrev[i]) ? devChange : dltDev;
            dltDev = ((0 < volt && volt <= devPrev[i]) || 6 * devPrev[i] < volt) ? -devChange : dltDev;
            IntVolt dev = devPrev[i] + dltDev;
            deviations[i] = (dev < minDev) ? minDev : dev; // clamp deviations at minDev
        }
    }

    void Detection::detection(const IntVolt *trace, const IntVolt *ref,
                              const IntVolt *baselines, const IntVolt *deviations,
                              IntChannel channelStart, IntChannel channelEnd, IntFrame t)
    {
        for (IntChannel i = channelStart; i < channelEnd; i++)
        {
            IntVolt volt = trace[i] - *ref - baselines[i]; // calc against updated baselines
            IntVolt dev = deviations[i];

            IntCalc voltThr = volt * thrQuant;
            IntCalc thr = threshold * dev;
            IntCalc minAvg = minAvgAmp * dev;
            IntCalc maxAHP = maxAHPAmp * dev;

            if (spikeTime[i] < 0) // not in spike
            {
                if (voltThr > thr) // threshold crossing
                {
                    spikeTime[i] = 0;
                    spikeAmp[i] = volt;
                    spikeArea[i] = voltThr;
                    hasAHP[i] = false;
                }
                continue;
            }
            // else: during a spike
            spikeTime[i]++;
            // 1 <= spikeTime[i]

            if (spikeTime[i] < ampAvgDur) // sum up area in ampAvgDur
            {
                spikeArea[i] += voltThr;
                if (spikeAmp[i] < volt) // larger amp found
                {
                    spikeTime[i] = 0; // reset peak to current
                    spikeAmp[i] = volt;
                    // but accumulate area (already added)
                    hasAHP[i] = false;
                }
                continue;
            }
            // else: ampAvgDur <= spikeTime[i]

            if (spikeTime[i] < spikeDur)
            {
                if (voltThr < maxAHP) // AHP found
                {
                    hasAHP[i] = true;
                }
                else if (spikeAmp[i] < volt) // larger amp found
                {
                    spikeTime[i] = 0; // reset peak to current
                    spikeAmp[i] = volt;
                    spikeArea[i] += voltThr; // but accumulate area
                    hasAHP[i] = false;
                }
                continue;
            }
            // else: spikeTime[i] == spikeDur, spike end

            if (spikeArea[i] > minAvg * ampAvgDur && // reach min area
                (hasAHP[i] || voltThr < maxAHP))     // AHP exist
            {
                pQueue->addSpike(Spike(t - spikeDur, i, spikeAmp[i]));
            }

            spikeTime[i] = -1; // reset counter even if not spike

        } // for i

    } // Detection::detection

} // namespace HSDetection
