#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <vector>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

namespace SpkDonline {
class Detection {
  int NChannels; // number of channels; is set when reading the data
  int *ChInd;    // indices for parallelization
  // Variables for variance and mean
  int *Qd; // noise amplitude
  int *Qm; // median
  // Variables for the spike detection
  int *Sl;      // counter for spike length
  bool *AHP;    // counter for repolarizing current
  int *Amp;     // buffers spike amplitude
  int *SpkArea; // integrates over spike
  // Parameters for variance and mean updates
  const int Tau_m0 = 4;  // timescale for updating Qm (increment is Qd/Tau_m)
  const int Qdmin = 200; // set minimum value of Qd
  // Parameters for spike detection
  int threshold; // threshold to detect spikes >11 is likely to be real
                     // spikes, but can and should be sorted afterwards
  int AHPthr; // signal should go below that threshold within MaxSl-Slmin frames
  int MaxSl;     // dead time in frames after peak, used for further testing
  int MinAvgAmp; // minimal avg. amplitude of peak (in units of Qd)
  int MinSl;     // length considered for determining avg. spike amplitude
  // Parameters for reading data
  long tInc; // 100, increment for reading data, has to be changed in main
             // program as well
  const int Ascale = -64; // factor to multiply to raw traces to increase
                          // resolution; definition of ADC counts had been
                          // changed!
  const int Voffset = 0;  // mean ADC counts, as initial value for Qm
  // Parameters for recalibration events and artefact handling
  const int artT = 10; // to use after artefacts; to update Qm for 10 frames
  int *A;              // control parameter for amplifier effects
  // Files to save the spikes etc.
  int Sampling;
  std::ofstream w; // for spikes
  std::ofstream wShapes; // for shapes
  int *Aglobal;
  int *Slice;
  int a; // buffer for Iterate()

  // Shapes print out parameters
  int** ChInd10;
  int Window = 31; // Change this to change size of cutout

public:
  Detection();
  ~Detection();
  void InitDetection(long nFrames, double nSec, int sf, int NCh, long ti,
                     long int *Indices, int agl);
  void SetInitialParams(int thres, int maa, int ahpthr, int maxsl, int minsl);
  void openSpikeFile(const char *name);
  void openFiles(const char *spikes, const char *shapes);
  void MedianVoltage(unsigned short *vm);
  void MeanVoltage(unsigned short *vm, int tInc, int tCut);
  void Iterate(unsigned short *vm, long t0, int tInc, int tCut, int tCut2);
  void FinishDetection();
};
};
