#ifndef SPIKEHANDLER_H
#define SPIKEHANDLER_H

#include "ProcessSpikes.h"

using namespace std;

// Define Infinite (Using INT_MAX caused overflow problems)
#define INF 10000

struct CustomLessThan
{
    bool operator()(tuple<int, float> const &lhs, tuple<int, float> const &rhs) const
    {
        return std::get<1>(lhs) < std::get<1>(rhs);
    }
};

namespace SpikeHandler {

void setInitialParameters(int _num_channels, int _spike_peak_duration, string file_name, \
						  int _noise_duration, float _noise_amp_percent, float _inner_radius, int* _masked_channels, float** _channel_positions, int** _neighbor_matrix, \
						  int _max_neighbors, int _num_com_centers, bool _to_localize, int _cutout_start, int _cutout_end, int _maxsl,
              bool _decay_filtering, bool _verbose);
void loadRawData(short *_raw_data, int _index_data, int _iterations, int maxFramesProcessed, int before_chunk, int after_chunk);
void setLocalizationParameters(int _aGlobal, int** _baselines, int _index_baselines);
void addSpike(int channel, int frame, int amplitude);
void terminateSpikeHandler();
//Inner neighbor creation methods
float channelsDist(int start_channel, int end_channel);
void fillNeighborLayerMatrices();
vector<int> getInnerNeighborsRadius(vector<tuple<int, float>> distances_neighbors, int central_channel);
int** createInnerNeighborMatrix();
int** createOuterNeighborMatrix();
Spike storeWaveformCutout(int cutout_size, Spike curr_spike);
Spike storeCOMWaveformsCounts(Spike curr_spike);

};

#endif
