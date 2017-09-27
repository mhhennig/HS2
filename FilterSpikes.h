#ifndef FILTERSPIKES_H  
#define FILTERSPIKES_H 

#include <deque>
#include "Parameters.h"

using namespace std;


namespace FilterSpikes {

bool areNeighbors(int channel_one, int channel_two);
Spike filterSpikes(Spike largest_amp_spike);

};

#endif