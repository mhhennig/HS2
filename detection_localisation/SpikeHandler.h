#ifndef SPIKEHANDLER_H
#define SPIKEHANDLER_H

#include "ProcessSpikes.h"

using namespace std;

// Define Infinite (Using INT_MAX caused overflow problems)
#define INF 10000

struct Point
{
    int x;
    int y;
    int channel;
};

struct Line
{
    Point p1;
    Point p2;
    float a; //ax
    float b; //by
    float c; //c
};

struct CustomLessThan
{
    bool operator()(tuple<int, float> const &lhs, tuple<int, float> const &rhs) const
    {
        return std::get<1>(lhs) < std::get<1>(rhs);
    }
};

namespace SpikeHandler {

void setInitialParameters(int _num_channels, int _spike_delay, int _spike_peak_duration, string file_name, \
						  int _noise_duration, float _noise_amp_percent, float _inner_radius, int* _masked_channels, int** _channel_positions, int** _neighbor_matrix, \
						  int _max_neighbors, bool _to_localize, int _cutout_start, int _cutout_end, int _maxsl, bool _verbose);
void loadRawData(short* _raw_data, int _index_data, int _iterations, int _frames, int _additional_data);
void setLocalizationParameters(int _aGlobal, int** _baselines, int _index_baselines);
void addSpike(int channel, int frame, int amplitude);
void terminateSpikeHandler();
//Inner neighbor creation methods
float channelsDist(int start_channel, int end_channel);
void fillNeighborLayerMatrices();
Line createLine(Point p1, Point p2);
void createBoundaryLines(Point curr_point, vector<Point> &boundary_points, vector<Line> &boundary_lines);
bool acceptAsBoundaryPoint(Point curr_point, Point central_point, vector<Line> &boundary_lines);
vector<int> getInnerNeighborsBounding(vector<tuple<int, float>> distances_neighbors, int central_channel);
vector<int> getInnerNeighborsRadius(vector<tuple<int, float>> distances_neighbors, int central_channel);
float distBetweenPoints(Point p1, Point p2);
int** createInnerNeighborMatrix();
int** createOuterNeighborMatrix();

};

#endif
