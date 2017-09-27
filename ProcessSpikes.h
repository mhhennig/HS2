#ifndef PROCESSSPIKES_H  
#define PROCESSSPIKES_H 

#include <iostream>
#include <stdio.h>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <deque>
#include "Parameters.h"
#include "FilterSpikes.h"


using namespace std;

namespace ProcessSpikes {

void filterSpikes(ofstream& spikes_filtered_file);
void filterLocalizeSpikes(ofstream& spikes_filtered_file);

};

#endif