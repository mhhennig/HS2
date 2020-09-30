import argparse

from herdingspikes.probe import NeuroPixel
from herdingspikes.hs2 import HSDetection, HSClustering

# raw data file
data_path = "rawDataSample.bin"

# detection parameters
file_directory = 'results/'
file_name = 'ProcessedSpikesNeuropixel'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Spike detection and clustering using Heardingspikes')
    parser.add_argument('-i', '--input', type=str,
        help="path to binary file to read from, default is rawDataSample.bin")

    args = parser.parse_args()

    if args.input:
        data_path = args.input
    print("Using " + data_path + "as data input")


    print("Initialisation ...")
    Probe = NeuroPixel(data_path, neighbor_radius=120,
                    masked_channels=[36, 75, 112, 151, 188, 227, 264, 303, 340, 379, 384])
    hs_detection = HSDetection(Probe, num_com_centers=2, threshold=20,
                    left_cutout_time=0.4, right_cutout_time=1.0,
                    maa=0, amp_evaluation_time=0.1, spk_evaluation_time=0.4,
                    ahpthr=0, out_file_name=file_name, 
                    file_directory_name=file_directory, decay_filtering=True, save_all=True)
    print("Detecting spikes ...")
    hs_detection.DetectFromRaw()
    hs_detection.LoadDetected()

    print("Clustering spikes ...")
    hs_clustering = HSClustering(hs_detection)
    hs_clustering.SaveHDF5(file_directory + file_name+"_sorted.hdf5")

