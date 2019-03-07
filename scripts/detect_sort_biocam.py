# Detect and sort multiple recordings
# this example reads Biocam files

import sys
import glob

# import locally (works without install)
sys.path.append('..')
sys.path.insert(0, '../herdingspikes/')
sys.path.insert(0, '../')

from herdingspikes.hs2 import HSDetection,HSClustering
from herdingspikes.probe import BioCam
from datetime import datetime
import getopt
import h5py
import sklearn

if __name__ == '__main__':
   
    sklearn.set_config(assume_finite=True, working_memory=1024*8)
 
    def printargs():
        print('HS2.py -i <inputfile or file mask> [-m <mbf> -a <0/1> -p <0/1>]')

    if len(sys.argv)==1:
        printargs()
        sys.exit(2)

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:", ["ifile="])
    except getopt.GetoptError:
        printargs()
        sys.exit(2)

    for opt, arg in opts:
        print(opt, arg)
        if opt == '-h':
            printargs()
            sys.exit()
        elif opt in ("-i", "--ifile"):
            filemask = arg

    try:
        sfs = glob.glob(filemask)
        print('File(s):')
        for f in sfs:
            print(f)
    except:
        print('cannot find ' + filemask)
        sys.exit(2)

    # detection parameters
    # these work for 7kHz Biocam/retina
    to_localize = True
    cutout_start = 4
    cutout_end = 20
    cutout_length = cutout_end+cutout_start+1
    threshold = 22
    num_com_centers = 2
    # what to do
    detect_it = False #True
    sort_it = True
            
    if sfs[0].find('.hdf5')>0:
        print('Re-clustering')
        detect_it = False
    else:
        print('Detecting and clustering')
            
    # Detect and localise
    if detect_it is True:
        print("detecting...")
        startTime = datetime.now()
        for f in sfs:
            out_file = f.replace('.brw', '')
            print("detecting "+f)
            # make a probe object first, and then set up the detector
            Probe = BioCam(f)
            H = HSDetection(Probe, to_localize, num_com_centers, cutout_start, cutout_end, threshold,
                        maa=0, maxsl=12, minsl=3, ahpthr=0, out_file_name=out_file,
                        save_all=False)
            H.DetectFromRaw()
            sampling = H.probe.fps
            # garbage collect
            del H
            del Probe
        print('Time taken for detection: ' + str(datetime.now() - startTime))
    else:
        # get sampling rate to store in clustered file
        print('getting sampling rate from '+sfs[5])
        f = h5py.File(sfs[5],'r')
        recVars = f.require_group('3BRecInfo/3BRecVars/')
        sampling = recVars['SamplingRate'].value[0]
        f.close()
        print('Sampling rate:'+str(sampling))

    # cluster everything
    if sort_it is True:
        print("reading files for clustering...")
        #if detect_it is True:
        out_files = [f.replace('.brw', '.bin') for f in sfs]
        #else:
        #    out_files = sfs
        C = HSClustering(out_files, cutout_length=cutout_length)
        C.ShapePCA(pca_ncomponents=2, pca_whiten=True)
        print("starting Mean Shift...")
        startTime = datetime.now()
        # parameters work for Biocam geometry, 42um pitch
        C.CombinedClustering(alpha=0.28, bandwidth=0.28, bin_seeding=False, cluster_subset=200000, n_jobs=-1)
        #if detect_it is True:
        sorted_files = [f.replace('.brw', '_HS2_clustered.hdf5') for f in sfs]
        #else:
        #    sorted_files = [f.replace('.hdf5', '_HS2_clustered.hdf5') for f in sfs]
        print('Time taken for sorting: ' + str(datetime.now() - startTime))

        # save all spikes in hdf5
        C.SaveHDF5(sorted_files,sampling=sampling)
