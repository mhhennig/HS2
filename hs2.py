from __future__ import division
import pandas as pd
import numpy as np
import h5py
from detection_localisation.detect import detectData
from matplotlib import pyplot as plt
#from sklearn.cluster import MeanShift # joblib things are broken
from clustering.mean_shift_ import MeanShift
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
    Ullo, S., ... & Hennig, M. H. (2017). Unsupervised spike sorting for
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
        self.shapecache = None
        self.IsClustered = False
        self.HasFeatures = False

    def Save_HDF5_legacy(self, filename, compression=None):
        """
        Saves data, cluster centres and ClusterIDs to a hdf5 file.
        Offers compression of the shapes, 'lzf'
        appears a good trade-off between speed and performance.
        """
        g = h5py.File(filename, 'w')
        g.create_dataset("data", data=np.vstack(
            (self.spikes.x, self.spikes.y)))
        g.create_dataset("Sampling", data=self.probe.fps)
        g.create_dataset("times", data=self.spikes.t)
        # g.create_dataset("expinds", data=self.__expinds)
        if self.IsClustered:
            g.create_dataset("centres", data=self.centerz.T)
            g.create_dataset("cluster_id", data=self.spikes.cl)
        # this is still a little slow (and perhaps memory intensive)
        # but I have not found a better way:
        sh_tmp = np.empty((self.cutout_length, self.spikes.Shape.size),
                          dtype=int)
        for i in range(self.spikes.Shape.size):
            sh_tmp[:, i] = self.spikes.Shape[i]
        g.create_dataset("shapes", data=sh_tmp, compression=compression)
        g.close()

    def Load_HDF5_legacy(self, filename, compute_amplitudes=False,
                         chunk_size=500000, compute_cluster_sizes=False):
        """
        Load data, cluster centres and ClusterIDs from a hdf5 file created with
        HS1.
        """
        g = h5py.File(filename, 'r')
        print('Reading from ' + filename)
        print('Creating memmapped cache for shapes, reading in chunks, ' +
              'converting to integer...')
        self.shapecache = np.memmap(
            "tmp.bin", dtype=np.int32, mode="w+", shape=g['shapes'].shape[::-1])
        for i in range(g['shapes'].shape[1] // chunk_size + 1):
            tmp = (1000*np.transpose(g['shapes'][:, i*chunk_size:(
                i+1)*chunk_size])).astype(np.int32)
            inds = np.where(tmp > 20000)[0]
            tmp[inds] = 0
            print('Found ' + str(len(inds)) +
                  ' data points out of linear regime in chunk '+str(i+1))
            self.shapecache[i * chunk_size:(i + 1) * chunk_size] = tmp[:]

        self.probe.fps = g['Sampling'].value
        self.centerz = g['centres'].value.T
        self.NClusters = len(self.centerz)
        self.cutout_length = self.shapecache.shape[1]

        self.spikes = pd.DataFrame({'ch': np.zeros(g['times'].shape[0],
                                                   dtype=int),
                                    't': g['times'],
                                    'Amplitude': np.zeros(g['times'].shape[0],
                                                          dtype=int),
                                    'x': g['data'][0, :],
                                    'y': g['data'][1, :],
                                    'Shape': list(self.shapecache),
                                    'cl': g['cluster_id']
                                    }, copy=False)

        print('Found ' + str(len(self.spikes.t)) +
              ' spikes in ' + str(self.NClusters) + ' units...')

        if compute_amplitudes:
            print('Computing amplitudes...')
            self.spikes.Amplitude = [np.min(s) for s in self.shapecache]
            _avgAmpl = [np.mean(self.spikes.Amplitude[cl])
                        for cl in range(self.NClusters)]
        else:
            _avgAmpl = np.zeros(self.NClusters, dtype=int)

        if compute_cluster_sizes:
            _cls = [np.sum(self.spikes.cl == cl)
                    for cl in range(self.NClusters)]
        else:
            _cls = np.zeros(self.NClusters, dtype=int)

        dic_cls = {'ctr_x': self.centerz[:, 0],
                   'ctr_y': self.centerz[:, 1],
                   'Color': 1. * np.random.permutation(
            self.NClusters) / self.NClusters,
            'Size': _cls,
            'AvgAmpl': _avgAmpl}
        self.clusters = pd.DataFrame(dic_cls)

        self.IsClustered = True
        self.HasFeatures = False
        g.close()
        print('...done.\nWarning: HS2 does not know raw data is not available.')

    def LoadDetected(self, file_name, cutout_length):
        """
        Reads a binary file with spikes detected with the DetectFromRaw() method
        """
        self.cutout_length = cutout_length
        sp_flat = np.memmap(file_name + ".bin", dtype=np.int32, mode="r")
        assert sp_flat.shape[0] // (cutout_length + 5) is not 1. * \
            sp_flat.shape[0] / (cutout_length + 5), \
            "spike data has wrong dimensions"
        self.shapecache = sp_flat.reshape((sp_flat.shape[0] // (
            cutout_length + 5), cutout_length + 5))
        self.spikes = pd.DataFrame({'ch': self.shapecache[:, 0],
                                    't': self.shapecache[:, 1],
                                    'Amplitude': self.shapecache[:, 2],
                                    'x': self.shapecache[:, 3] / 1000,
                                    'y': self.shapecache[:, 4] / 1000,
                                    'Shape': list(self.shapecache[:, 5:])
                                    }, copy=False)
        self.IsClustered = False
        self.HasFeatures = False
        print('Detected and read ' + str(self.spikes.shape[0]) + ' spikes.')

    def DetectFromRaw(self, to_localize, cutout_start, cutout_end, threshold,
                      maa=0, maxsl=12, minsl=3, ahpthr=0, tpre=1.0, tpost=2.2,
                      out_file_name="ProcessedSpikes.bin", load=True):
        """
        This function is a wrapper of the C function `detectData`. It takes
        the raw data file, performs detection and localisation, saves the result
        to `ProcessedSpikes` and loads the latter into memory by calling
        `LoadDetected`.

        Arguments:
        file_name -- the path to the raw data file.
        to_localize
        cutout_start
        cutout_end
        threshold
        maa
        maxsl
        minsl
        ahpthr
        out_file_name
        load -- load the detected spikes when finished
        """
        detectData(self.probe, str.encode(out_file_name),
                   to_localize, self.probe.fps, threshold,
                   cutout_start, cutout_end,
                   maa, maxsl, minsl, ahpthr, tpre, tpost)
        if load:
            # reload data into memory
            cutout_length = cutout_start + cutout_end + 1
            self.LoadDetected(out_file_name, cutout_length)

    def ShapePCA(self, pca_ncomponents=2, pca_whiten=True, chunk_size=1000000):
        pca = PCA(n_components=pca_ncomponents, whiten=pca_whiten)
        if self.spikes.shape[0] > 1e6:
            print("Fitting PCA using 1e6 out of",
                  self.spikes.shape[0], "spikes...")
            inds = np.random.choice(self.spikes.shape[0], int(1e6),
                                    replace=False)
            pca.fit(np.array(list(self.spikes.Shape[inds])))
        else:
            print("Fitting PCA using " + str(self.spikes.shape[0]) + " spikes...")
            pca.fit(np.array(list(self.spikes.Shape)))
        self.pca = pca
        self.HasFeatures = True
        _pcs = np.empty((self.spikes.shape[0], pca_ncomponents))
        for i in range(self.spikes.shape[0] // chunk_size + 1):
            _pcs[i*chunk_size:(i + 1)*chunk_size, :] = pca.transform(np.array(
                list(self.spikes.Shape[i * chunk_size:(i + 1) * chunk_size])))

        return _pcs

    def CombinedClustering(self, alpha, clustering_algorithm=MeanShift,
                           pca_ncomponents=2, pca_whiten=True,
                           recompute_pca=False, **kwargs):
        """
        Performs PCA on the available spike shapes, and clusters spikes based
        on their (x, y) location and on the principal components of the shape.
        Cluster memberships are available as self.spikes.cl.
        Cluster information is available in the self.clusters dataframe.

        Arguments:
        alpha -- the weight given to PCA components, relative to
        spatial components
        clustering_algorithm -- a sklearn.cluster class, defaults to
        sklearn.cluster.MeanShift. sklearn.cluster.DBSCAN was also tested.
        pca_ncomponents -- number of PCA components to be considered (def. 2).
        pca_whiten -- whiten option of the PCA algorithm.
        **kwargs -- additional arguments are passed to the clustering class.
        This may include n_jobs > 1 for parallelisation.
        """

        if (not self.HasFeatures) or recompute_pca:
            self.cutouts_pca = self.ShapePCA(pca_ncomponents, pca_whiten)

        fourvec = np.vstack(([self.spikes.x], [self.spikes.y],
                             alpha * self.cutouts_pca.T)).T

        print('Clustering...')
        clusterer = clustering_algorithm(**kwargs)
        clusterer.fit(fourvec)
        self.spikes.cl = clusterer.labels_
        self.NClusters = len(np.unique(clusterer.labels_))
        # assert np.max(self.spikes.cl)+1 == self.NClusters
        print("Number of estimated clusters:", self.NClusters)
        in_cl = []
        for i in np.arange(self.NClusters):
            in_cl.append(list(np.where(self.spikes.cl == i)[0]))
        self.fourvec = fourvec
        centers = np.asarray([np.mean(fourvec[cl], axis=0) for cl in in_cl])
        self.centerz = centers
        self.in_cl = in_cl
        dic_cls = {'ctr_x': centers[:, 0],
                   'ctr_y': centers[:, 1],
                   'Color': 1. * np.random.permutation(
            self.NClusters) / self.NClusters,
            'Size': [len(cl) for cl in in_cl],
            'AvgAmpl': [np.mean(self.spikes.Amplitude[cl])
                        for cl in in_cl]}
        self.clusters = pd.DataFrame(dic_cls)
        self.IsClustered = True

    def PlotTracesChannels(self, eventid, ax=None,
                           window_size=200, cutout_start=6):
        """
        Draw a figure with an electrode and its neighbours, showing the raw
        traces and events. Note that this requires loading the raw data in
        memory again.

        Arguments:
        eventid -- centers, spatially and temporally, the plot to a specific
        event id.
        ax -- a matplotlib axes object where to draw. Defaults to current axis.
        window_size -- number of samples shown around a spike
        cutout_start -- n. frames recorded before the spike peak in the cutout
        """
        pos, neighs = self.probe.positions, self.probe.neighbors

        event = self.spikes.loc[eventid]
        print("Spike detected at channel: ", event.ch)
        print("Spike detected at frame: ", event.t)
        cutlen = len(event.Shape)
        assert window_size > cutlen, "window_size is too small"
        dst = np.abs(pos[event.ch][0] - pos[neighs[event.ch]][:, 0])
        interdistance = np.min(dst[dst > 0])
        if ax is None:
            ax = plt.gca()

        # scatter of the large grey balls for electrode location
        plt.scatter(np.array(pos)[neighs[event.ch], 0],
                    np.array(pos)[neighs[event.ch], 1],
                    s=1600, alpha=0.2)
        x = pos[[neighs[event.ch], 0]]
        y = pos[[neighs[event.ch], 1]]
        for i, txt in enumerate(neighs[event.ch]):
            ax.annotate(txt, (x[i] + 1, y[i] + 3))

        ws = window_size // 2
        t1 = np.max((0, event.t - ws))
        t2 = event.t + ws
        scale = interdistance/110.
        trange = (np.arange(t1, t2)-event.t)*scale
        start_bluered = event.t-t1-cutout_start
        trange_bluered = trange[start_bluered:start_bluered+cutlen]
        trange_bluered = np.arange(-cutout_start, -cutout_start+cutlen)*scale

        data = self.probe.Read(t1, t2).reshape((t2-t1, self.probe.num_channels))
        for n in neighs[event.ch]:
            if n not in self.probe.masked_channels:
                plt.plot(pos[n][0] + trange,
                         pos[n][1] + data[:, n]*scale, 'gray')
                plt.plot(pos[n][0] + trange_bluered,
                         pos[n][1] + data[start_bluered:start_bluered+cutlen,
                         n]*scale, 'b')
            else:
                plt.plot(pos[n][0] + trange,
                         pos[n][1] + data[:, n]*scale, 'gray')
                plt.plot(pos[n][0] + trange_bluered,
                         pos[n][1] + data[start_bluered:start_bluered+cutlen,
                         n]*scale, 'g')

        plt.plot(pos[event.ch][0] + trange_bluered,
                 pos[event.ch][1] + event.Shape * scale, 'r')
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

    def PlotAll(self, invert=False, show_labels=False, ax=None,
                max_show=200000, **kwargs):
        """
        Plots all the spikes currently stored in the class, in (x, y) space.
        If clustering has been performed, each spike is coloured according to
        the cluster it belongs to.

        Arguments:
        invert -- (boolean, optional) if True, flips x and y
        show_labels -- (boolean, optional) if True, annotates each cluster
        centre with its cluster ID.
        ax -- a matplotlib axes object where to draw. Defaults to current axis.
        max_show -- maximum number of spikes to show
        **kwargs -- additional arguments are passed to pyplot.scatter
        """
        if ax is None:
            ax = plt.gca()
        x, y = self.spikes.x, self.spikes.y
        if invert:
            x, y = y, x
        if self.spikes.shape[0] > max_show:
            inds = np.random.choice(
                self.spikes.shape[0], max_show, replace=False)
            print(
                'We have ' + str(self.spikes.shape[0]) +
                ' spikes, only showing ' + str(max_show))
        else:
            inds = np.arange(self.spikes.shape[0])
        c = plt.cm.hsv(self.clusters.Color[self.spikes.cl]) \
            if self.IsClustered else 'r'
        ax.scatter(x[inds], y[inds], c=c[inds], **kwargs)
        if show_labels and self.IsClustered:
            ctr_x, ctr_y = self.clusters.ctr_x, self.clusters.ctr_y
            if invert:
                ctr_x, ctr_y = ctr_y, ctr_x
            for cl in range(self.NClusters):
                if ~np.isnan(ctr_y[cl]):  # hack, why NaN positions in DBScan?
                    ax.annotate(str(cl), [ctr_x[cl], ctr_y[cl]], fontsize=16)
                    # seems this is a problem when zooming with x/ylim
        return ax

    def PlotShapes(self, units, nshapes=100, ncols=4, ax=None, ylim=None):
        """
        Plot a sample of the spike shapes contained in a given set of clusters
        and their average.

        Arguments:
        units -- a list of the cluster IDs to be considered.
        nshapes -- the number of shapes to plot (default 100).
        ncols -- the number of columns under which to distribute the plots.
        ax -- a matplotlib axis object (defaults to current axis).
        """
        nrows = np.ceil(len(units) / ncols)
        if ax is None:
            plt.figure(figsize=(3 * ncols, 3 * nrows))
        cutouts = np.array(list(self.spikes.Shape))

        # all this is to determine suitable ylims
        if ylim is None:
            meanshape = np.mean(cutouts, axis=0)
            maxy, miny = meanshape.max(), meanshape.min()
            varshape = np.var(cutouts, axis=0)
            varmin = varshape[np.argmin(meanshape)]
            varmax = varshape[np.argmax(meanshape)]
            maxy += 4. * np.sqrt(varmax)
            miny -= 2. * np.sqrt(varmin)
            ylim = [miny, maxy]

        for i, cl in enumerate(units):
            inds = np.where(self.spikes.cl == cl)[0]
            if ax is None:
                plt.subplot(nrows, ncols, i + 1)
            plt.plot(cutouts[inds[:100], :].T, 'gray')
            plt.plot(np.mean(cutouts[inds, :], axis=0),
                     c=plt.cm.hsv(self.clusters.Color[cl]), lw=4)
            plt.ylim(ylim)
            plt.title("Cluster " + str(cl))
