#ifndef FILTERSPIKES_H
#define FILTERSPIKES_H

#include "Parameters.h"

using namespace std;


namespace FilterSpikes {

Spike filterSpikes(Spike largest_amp_spike, ofstream& filteredsp);
Spike findMaxSpikeNeighbor(Spike first_spike);
Spike updateNeighborsMaxSpike(Spike max_spike);
void filterOuterNeighbors(Spike max_spike, ofstream& filteredsp);
bool filteredOuterSpike(Spike outer_spike, Spike max_spike);
int getClosestInnerNeighborChannel(int outer_channel, int central_channel);
Spike getSpikefromChannel(int spike_channel);
Event getEventfromChannel(int channel, Spike curr_spike, bool is_inner_event);
void filterInnerNeighbors(Spike max_spike, ofstream& filteredsp);
float channelsDist(int start_channel, int end_channel);
bool isInnerNeighbor(int central_channel, int curr_channel);
bool areNeighbors(int channel_one, int channel_two);
float posToNegRatio(Spike spike);
float repolarizationTime(Spike spike);
double areaUnderSpike(Spike spike);

};

#endif
