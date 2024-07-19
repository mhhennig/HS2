#include <algorithm>

#include "SpikeQueue.h"
#include "Detection.h"
#include "QueueProcessor/FirstElemProcessor.h"
#include "QueueProcessor/MaxSpikeFinder.h"
#include "QueueProcessor/SpikeDecayFilterer.h"
#include "QueueProcessor/SpikeFilterer.h"
#include "SpikeProcessor/SpikeLocalizer.h"
#include "SpikeProcessor/SpikeShapeWriter.h"

using namespace std;

#define pushFirstElemProc(pSpkProc)                  \
    do                                               \
    {                                                \
        spkProcs.push_back(pSpkProc);                \
        pQueProc = new FirstElemProcessor(pSpkProc); \
        queProcs.push_back(pQueProc);                \
    } while (false)

namespace HSDetection
{
    SpikeQueue::SpikeQueue(Detection *pDet)
        : spikes((Spike *)new char[pDet->chunkSize * pDet->numChannels * sizeof(Spike)]), spikeCnt(0),
          queue(), queProcs(), spkProcs(), pRresult(&pDet->result),
          procDelay(max(pDet->cutoutEnd - pDet->spikeDur, pDet->riseDur) + pDet->temporalJitter + 1)
    {
        SpikeProcessor *pSpkProc;
        QueueProcessor *pQueProc;

        pQueProc = new MaxSpikeFinder(&pDet->probeLayout, pDet->temporalJitter);
        queProcs.push_back(pQueProc);

        if (pDet->decayFilter)
        {
            pQueProc = new SpikeDecayFilterer(&pDet->probeLayout, pDet->temporalJitter, pDet->decayRatio);
        }
        else
        {
            pQueProc = new SpikeFilterer(&pDet->probeLayout, pDet->temporalJitter);
        }
        queProcs.push_back(pQueProc);

        if (pDet->localize)
        {
            pSpkProc = new SpikeLocalizer(&pDet->probeLayout, &pDet->trace, &pDet->commonRef, &pDet->runningBaseline,
                                          pDet->temporalJitter, pDet->riseDur);
            pushFirstElemProc(pSpkProc);
        }

        if (pDet->saveShape)
        {
            pSpkProc = new SpikeShapeWriter(pDet->filename, &pDet->trace, pDet->cutoutStart, pDet->cutoutEnd);
            pushFirstElemProc(pSpkProc);
        }
    }

    SpikeQueue::~SpikeQueue()
    {
        // should release QueueProc first because SpikeProc can be wrapped inside
        for_each(queProcs.begin(), queProcs.end(),
                 [](QueueProcessor *pQueProc)
                 { delete pQueProc; });
        for_each(spkProcs.begin(), spkProcs.end(),
                 [](SpikeProcessor *pSpkProc)
                 { delete pSpkProc; });

        delete[](char *) spikes;
    }

    void SpikeQueue::procFront()
    {
        for_each(queProcs.begin(), queProcs.end(),
                 [this](QueueProcessor *pQueProc)
                 { (*pQueProc)(this); });

        pRresult->push_back(std::move(*queue.begin()));
        queue.erase(queue.begin());
    }

    void SpikeQueue::process()
    {
        sort(spikes, spikes + spikeCnt,
             [](const Spike &lhs, const Spike &rhs)
             { return lhs.frame < rhs.frame || (lhs.frame == rhs.frame && lhs.channel < rhs.channel); });

        for (IntResult i = 0; i < spikeCnt; i++)
        {
            while (!queue.empty() && queue.front().frame < spikes[i].frame - procDelay)
            {
                procFront();
            }

            queue.push_back(std::move(spikes[i]));
        }

        spikeCnt = 0; // reset for next chunk
    }

    void SpikeQueue::finalize()
    {
        while (!queue.empty())
        {
            procFront();
        }
    }

} // namespace HSDetection
