#ifndef LOCALIZESPIKES_H  
#define LOCALIZESPIKES_H

#include "Parameters.h"

using namespace std;


namespace LocalizeSpikes {

tuple<float,float> centerOfMass(vector<float>* amps, int spike_channel);
tuple<float,float> localizeSpike(Spike spike_to_be_localized);

};

#endif