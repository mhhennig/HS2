import numpy as np
import pandas as pd
from detect import detectData
from matplotlib import pyplot as plt
from sklearn.cluster import MeanShift
from sklearn.decomposition import PCA


class herdingspikes(object):
    """
    This class provides a simple interface to the detection, localisation and
    clustering of spike data from dense multielectrode arrays according to the
    methods described in the following papers:

    Muthmann, J. O., Amin, H., Sernagor, E., Maccione, A., Panas, D.,
    Berdondini, L., ... & Hennig, M. H. (2015). Spike detection for large neural
    populations using high density multielectrode arrays.
    Frontiers in neuroinformatics, 9.

    Hilgen, G., Sorbaro, M., Pirmoradian, S., Muthmann, J. O., Kepiro, I. E.,
    Ullo, S., ... & Murino, V. (2017). Unsupervised spike sorting for
    large-scale, high-density multielectrode arrays.
    Cell reports, 18(10), 2521-2532.

    Usage:
    1. call the constructor with a single argument -- a NeuralProbe object:
        HS = herdingspikes(Probe)
    2. You can now load your data in three ways: with LoadDetected()
    (automatically loads ProcessedSpikes); with DetectFromRaw (takes a path
    to the raw data and all detection parameters) or with LoadH5 (takes a
    previously saved instance of this class)
    """
    def __init__(self, probe):
        self.probe = probe

    # def SaveH5(self, filename):
    #     store = pd.HDFStore(filename)
    #     store['spikes'] = self.spikes
    #     if self.IsClustered:
    #         store['clusters'] = self.clusters
    #     store.close()
    #
    # def LoadH5(self, filename):
    #     store = pd.HDFStore(filename)
    #     self.spikes = store['spikes']
    #     if 'clusters' in store.keys():
    #         self.clusters = store['clusters']
    #         self.IsClustered = True
    #         self.NClusters = len(self.clusters)
    #         assert np.max(self.spikes.cl)+1 == self.NClusters
    #     else:
    #         self.IsClustered = False
    #     store.close()

    def LoadDetected(self):
        """
        Reads the `ProcessedSpikes` file present in the current directory.
        """
        sp = np.loadtxt('ProcessedSpikes')
        self.spikes = pd.DataFrame({'ch': sp[:, 0].astype(int),
                                    't': sp[:, 1].astype(int),
                                    'Amplitude': sp[:, 2],
                                    'x': sp[:, 3],
                                    'y': sp[:, 4],
                                    'Shape': list(sp[:, 5:])
                                    })
        self.IsClustered = False

    def DetectFromRaw(self, datapath,
                      to_localize, cutout_start, cutout_end, threshold,
                      maa=0, maxsl=12, minsl=3, ahpthr=0, data_format='flat'):
        """
        This function is a wrapper of the C function `detectData`. It takes
        the raw data file, performs detection and localisation, saves the result
        to `ProcessedSpikes` and loads the latter into memory by calling
        `LoadDetected`.

        Arguments:
        datapath -- the path to the raw data file.
        to_localize
        cutout_start
        cutout_end
        threshold
        maa
        maxsl
        minsl
        ahpthr
        """
        probe = self.probe
        detectData(datapath, probe.positions_file_path,
                   probe.neighbors_file_path,
                   probe.num_channels, probe.num_recording_channels,
                   probe.spike_delay, probe.spike_peak_duration,
                   probe.noise_duration, probe.noise_amp_percent,
                   probe.max_neighbors,
                   to_localize, probe.fps, threshold, cutout_start, cutout_end,
                   maa, maxsl, minsl, ahpthr, data_format=probe.data_format)
        # reload data into memory
        self.LoadDetected()

    def PlotTracesChannels(self, datapath, eventid, ax=None):
        """
        Draw a figure with an electrode and its neighbours, showing the raw
        traces and events. Note that this requires loading the raw data in
        memory again.

        Arguments:
        datapath -- the path to the raw data file.
        eventid -- centers, spatially and temporally, the plot to a specific
        event id.
        ax -- a matplotlib axes object where to draw. Defaults to current axis.
        """
        if ax is None:
            ax = plt.gca()
        pos, neighs = self.probe.positions, self.probe.neighbors
        data = np.fromfile(datapath,
                           dtype=np.int16).reshape((-1,
                                                    self.probe.num_channels))

        event = self.spikes.loc[eventid]
        print("Spike detected at channel: ", event.ch)
        print("Spike detected at frame: ", event.t)
        cutlen = len(event.Shape)

        plt.scatter(np.array(pos)[neighs[event.ch], 0],
                    np.array(pos)[neighs[event.ch], 1],
                    s=1600, alpha=0.2)

        t1 = np.max((0, event.t - 100))
        t2 = event.t + 100

        for n in neighs[event.ch]:
            plt.plot(pos[n][0]+np.arange(t2-t1)/10.-10,
                     pos[n][1] + data[t1:t2, n]/10., 'gray')
            plt.plot(pos[n][0]+(np.arange(cutlen)+event.t-t1-10)/10.-10,
                     pos[n][1] + data[event.t-10:event.t+cutlen-10, n]/10.,
                     'b')
        plt.plot(pos[event.ch][0]+(np.arange(cutlen)+event.t-t1-10)/10.-10,
                 pos[event.ch][1] + event.Shape/10., 'r')
        plt.scatter(event.x, event.y, s=80, c='r')

    def PlotDensity(self, binsize=1., invert=False, ax=None):
        raise NotImplementedError()
        if ax is None:
            ax = plt.gca()
        x, y = self.spikes.x, self.spikes.y
        if invert:
            x, y = y, x
        binsx = np.arange(x.min(), x.max(), binsize)
        binsy = np.arange(y.min(), y.max(), binsize)
        h, xb, yb = np.histogram2d(x, y, bins=[binsx, binsy])
        ax.imshow(np.log10(h), extent=[xb.min(), xb.max(), yb.min(), yb.max()],
                  interpolation='none', origin='lower')
        return h, xb, yb

    def PlotAll(self, invert=False, show_labels=False, ax=None, **kwargs):
        """
        Plots all the spikes currently stored in the class, in (x, y) space.
        If clustering has been performed, each spike is coloured according to
        the cluster it belongs to.

        Arguments:
        invert -- (boolean, optional) if True, flips x and y
        show_labels -- (boolean, optional) if True, annotates each cluster
        centre with its cluster ID.
        ax -- a matplotlib axes object where to draw. Defaults to current axis.
        **kwargs -- additional arguments are passed to pyplot.scatter
        """
        if ax is None:
            ax = plt.gca()
        x, y = self.spikes.x, self.spikes.y
        if invert:
            x, y = y, x
        c = plt.cm.hsv(self.clusters.Color[self.spikes.cl]) \
            if self.IsClustered else 'r'
        ax.scatter(x, y, c=c, **kwargs)
        if show_labels and self.IsClustered:
            ctr_x, ctr_y = self.clusters.ctr_x, self.clusters.ctr_y
            if invert:
                ctr_x, ctr_y = ctr_y, ctr_x
            for cl in range(self.NClusters):
                ax.annotate(str(cl), [ctr_x[cl], ctr_y[cl]], fontsize=16)

    def CombinedClustering(self, alpha, clustering_algorithm=MeanShift,
                           pca_ncomponents=2, pca_whiten=True,
                           **kwargs):
        """
        Performs PCA on the available spike shapes, and clusters spikes based
        on their (x, y) location and on the principal components of the shape.
        Cluster memberships are available as self.spikes.cl.
        Cluster information is available in the self.clusters dataframe.

        Arguments:
        alpha -- the weight given to PCA components (spatial components are
        assigned weight 1).
        clustering_algorithm -- a sklearn.cluster class, defaults to
        sklearn.cluster.MeanShift. sklearn.cluster.DBSCAN was also tested.
        pca_ncomponents -- number of PCA components to be considered (def. 2).
        pca_whiten -- whiten option of the PCA algorithm.
        **kwargs -- additional arguments are passed to the clustering class.
        This may include n_jobs > 1 for parallelisation.
        """
        pca = PCA(n_components=pca_ncomponents, whiten=pca_whiten)
        cutouts_pca = pca.fit_transform(np.array(list(self.spikes.Shape)))
        fourvec = np.vstack(([self.spikes.x], [self.spikes.y],
                             alpha*cutouts_pca.T)).T

        clusterer = clustering_algorithm(**kwargs)
        clusterer.fit(fourvec)
        self.spikes.cl = clusterer.labels_
        self.NClusters = len(np.unique(clusterer.labels_))
        # assert np.max(self.spikes.cl)+1 == self.NClusters
        print("Number of estimated clusters:", self.NClusters)

        in_cl = []
        for i in np.arange(self.NClusters):
            in_cl.append(list(np.where(self.spikes.cl == i)[0]))
        centers = np.asarray([np.mean(fourvec[cl], axis=0) for cl in in_cl])
        dic_cls = {'ctr_x': centers[:, 0],
                   'ctr_y': centers[:, 1],
                   'Color': 1.*np.random.permutation(
                        self.NClusters)/self.NClusters,
                   'Size': [len(cl) for cl in in_cl],
                   'AvgAmpl': [np.mean(self.spikes.Amplitude[cl])
                               for cl in in_cl]}
        self.clusters = pd.DataFrame(dic_cls)
        self.IsClustered = True

    def PlotShapes(self, units, nshapes=100, ncols=4):
        """
        Plot a sample of the spike shapes contained in a given set of clusters
        and their average.

        Arguments:
        units -- a list of the cluster IDs to be considered.
        nshapes -- the number of shapes to plot (default 100).
        ncols -- the number of columns under which to distribute the plots.
        """
        nrows = np.ceil(len(units)/ncols)
        plt.figure(figsize=(3*ncols, 3*nrows))
        cutouts = np.array(list(self.spikes.Shape))
        for i, cl in enumerate(units):
            inds = np.where(self.spikes.cl == cl)[0]
            plt.subplot(nrows, ncols, i+1)
            plt.plot(cutouts[inds[:100], :].T, 'gray')
            plt.plot(np.mean(cutouts[inds, :], axis=0),
                     c=plt.cm.hsv(self.clusters.Color[cl]), lw=4)
            plt.ylim((-130, 90))
            plt.title("Cluster "+str(cl))
