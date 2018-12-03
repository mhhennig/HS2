#ifndef LOCALIZESPIKES_H
#define LOCALIZESPIKES_H

#include "Parameters.h"

using namespace std;


namespace LocalizeSpikes {

tuple<float, float> centerOfMass(deque<tuple<int, int>> centered_amps);
tuple<float, float> localizeSpike(Spike spike_to_be_localized);
tuple<float, float> reweightedCenterOfMass(deque<tuple<tuple<float, float>, int>> com_positions_amps);

};

#endif
