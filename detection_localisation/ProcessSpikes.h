#ifndef PROCESSSPIKES_H  
#define PROCESSSPIKES_H

#include "Parameters.h"
#include "FilterSpikes.h"
#include "LocalizeSpikes.h"

using namespace std;

namespace ProcessSpikes {

void filterSpikes(ofstream& spikes_filtered_file);
void filterLocalizeSpikes(ofstream& spikes_filtered_file);

};

#endif