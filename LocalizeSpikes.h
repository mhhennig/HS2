#ifndef LOCALIZESPIKES_H  
#define LOCALIZESPIKES_H

#include "Parameters.h"

using namespace std;


namespace LocalizeSpikes {

tuple<int, int> centerOfMass(int spike_channel);
tuple<int, int> localizeSpike(Spike spike_to_be_localized, int baseline_frame);

};

#endif