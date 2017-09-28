#ifndef PROCESSSPIKES_H  
#define PROCESSSPIKES_H

#include "Parameters.h"
#include "FilterSpikes.h"
#include "LocalizeSpikes.h"
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */


using namespace std;

namespace ProcessSpikes {

void filterSpikes(ofstream& spikes_filtered_file);
void filterLocalizeSpikes(ofstream& spikes_filtered_file);

};

#endif