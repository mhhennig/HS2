from hs2 import herdingspikes
from probe import NeuroPixel
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from random import shuffle
import random
from sklearn import svm
from sklearn import preprocessing
import numpy as np
import heapq

def getClosestClusters(cluster_id, num_neighbors, H):
    cluster_distances = [(cl_id, np.sqrt(((center - H.centerz[cluster_id])**2).sum())) for cl_id, center in enumerate(H.centerz)]
    closest_clusters = heapq.nsmallest(num_neighbors,cluster_distances, key=lambda X: X[1])
    return closest_clusters

def getRepresentativeWaveforms(cluster_id, num_waveforms, H):
    cluster = H.in_cl[cluster_id]
    spikes_dist_from_cluster = [(spike, np.sqrt(((H.fourvec[spike] - H.centerz[cluster_id])**2).sum())) for spike in cluster]
    representative_spikes = heapq.nsmallest(num_waveforms, spikes_dist_from_cluster, key=lambda X: X[1])
    representative_cutouts = np.array(list(H.spikes.Shape[[spike[0] for spike in representative_spikes]]))
    return (cluster_id, representative_cutouts, representative_spikes)

def createTrainingSet(cluster_representatives, H, pca):
    X = []
    Y = []
    for cluster in cluster_representatives:
        cl_id = cluster[0]
        cutouts = cluster[1]
        for i, example in enumerate(pca.transform(np.array(list(cutouts)))):
            spike_id = cluster[2][i][0]
            np.append(example, [H.spikes.x[spike_id]])
            np.append(example, [H.spikes.y[spike_id]])
            X.append(example)
            Y.append(cl_id)

    c = list(zip(X, Y))
    random.shuffle(c)
    X, Y = zip(*c)
    return (X, Y)

def getAllSpikeShapesinCluster(cluster_id, H):
    return np.array(list([H.spikes.Shape[spike] for spike in H.in_cl[cluster_id]]))

def getAllSpikePositionsinCluster(cluster_id, H):
    return [(H.spikes.x[spike], H.spikes.y[spike]) for spike in H.in_cl[cluster_id]]

def getProbabilitiesinClusterSVM(testing_set, testing_positions, labels, clf, pca):
    correct = 0
    wrong = 0
    probabilities = []
    testing_set_transforms = pca.transform(testing_set)
    for i, transform in enumerate(testing_set_transforms):
        np.append(transform, [testing_positions[i][0]])
        np.append(transform, [testing_positions[i][1]])
    for label in labels:
        for prediction in clf.predict(testing_set_transforms):
            if(prediction == label):
                correct += 1
            if(prediction != label):
                wrong += 1
        probabilities.append((label, correct/(correct + wrong)))
        correct = 0
        wrong = 0
    return probabilities

def main():
    data_path = '/disk/scratch/Cole/Hopkins_20160722_g0_t0.imec.ap_CAR.bin'
    to_localize = True
    cutout_start = 10
    cutout_end = 30
    threshold = 12
    masking = None
    file_name = 'ProcessedSpikes'

    Probe = NeuroPixel(data_file_path=data_path, fps=30000, masked_channels=masking)
    H = herdingspikes(Probe)

    H.DetectFromRaw(to_localize, cutout_start, cutout_end, threshold,
                    maa=0, maxsl=12, minsl=3, ahpthr=0)

    H.CombinedClustering(alpha=40,
                        bandwidth=20, bin_seeding=True, min_bin_freq=150,
                        pca_ncomponents=2, pca_whiten=True,
                        n_jobs=-1)

    H.Save_HDF5_legacy('out.hdf5')


    pca_whiten = True
    pca = PCA(n_components=10, whiten=pca_whiten)
    pca.fit(np.array(list(H.spikes.Shape)))

    cluster_representatives = []

    file = open('ClusterProbabilities.txt', 'w')

    #Get closest clusters to 0 (Minowski Distance)
    for i in range(20):
        closest_clusters = getClosestClusters(i, 5, H)

        #Get all representative waveforms from all neighbors
        for cluster in closest_clusters:
            representative_waveforms = getRepresentativeWaveforms(cluster[0], 1000, H)
            cluster_representatives.append(representative_waveforms)


        #Getting Training Data from all neighbors
        X, Y = createTrainingSet(cluster_representatives, H, pca)

        #Create and Train the classifier
        clf = svm.SVC(kernel = 'rbf')
        clf.fit(X , Y)

        #Run cluster 0 through the classifier and return probabilities
        testing_set = getAllSpikeShapesinCluster(i, H)
        positions_set = getAllSpikePositionsinCluster(i, H)
        labels = [cluster_id for cluster_id, _ in closest_clusters]
        probabilities = getProbabilitiesinClusterSVM(testing_set, positions_set, labels, clf, pca)
        probabilities.sort(key=lambda X: -X[1])
        cluster_representatives = []

        file.write("Cluster " + str(i))
        for probability in probabilities:
            file.write("%s\n" % probability[0])
            file.write("%s\n" % probability[1])


if __name__ == "__main__":
    main()
