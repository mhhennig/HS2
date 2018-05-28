#ifndef PROCESSSPIKES_H
#define PROCESSSPIKES_H

#include "Parameters.h"
#include "FilterSpikes.h"
#include "LocalizeSpikes.h"

using namespace std;

namespace ProcessSpikes {

void filterSpikes(ofstream& spikes_filtered_file, ofstream& filteredsp);
void filterLocalizeSpikes(ofstream& spikes_filtered_file, ofstream& filteredsp);
};

#endif
