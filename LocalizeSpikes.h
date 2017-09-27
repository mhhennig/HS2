#ifndef LOCALIZESPIKES_H  
#define LOCALIZESPIKES_H

#include <deque>
#include <tuple>
#include "Parameters.h"

using namespace std;


namespace FilterSpikes {

tuple<float,float> centerOfMass(tuple<int,int> channel_amps, int numNeighbors);
tuple<float,float> localizeSpike(Spike spike_to_be_localized);

};

#endif