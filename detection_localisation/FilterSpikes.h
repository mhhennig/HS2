#ifndef FILTERSPIKES_H  
#define FILTERSPIKES_H 

#include "Parameters.h"

using namespace std;


namespace FilterSpikes {

bool areNeighbors(int channel_one, int channel_two);
void eliminateDuplicates(Spike max_spike);
Spike filterSpikes(Spike largest_amp_spike);

};

#endif