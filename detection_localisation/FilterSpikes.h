#ifndef FILTERSPIKES_H
#define FILTERSPIKES_H

#include "Parameters.h"

using namespace std;


namespace FilterSpikes {

Spike filterSpikes(Spike largest_amp_spike, ofstream& filteredsp);
Spike findMaxSpikeNeighbor(Spike first_spike);
Spike updateNeighborsMaxSpike(Spike max_spike);
Spike updateOuterNeighbors(Spike max_spike, vector<tuple<int, int>> outer_spiking_neighbors);
Spike updateInnerNeighbors(Spike max_spike, vector<tuple<int, int>> inner_spiking_neighbors);
void filterOuterNeighbors(Spike max_spike, ofstream& filteredsp);
bool filteredOuterSpike(Spike outer_spike, Spike max_spike);
void filterInnerNeighbors(Spike max_spike, ofstream& filteredsp);
float channelsDist(int start_channel, int end_channel);
bool areNeighbors(int channel_one, int channel_two);

};

#endif
