import numpy as np
import pandas as pd
from detect import detectData
from matplotlib import pyplot as plt
from sklearn.cluster import MeanShift
from sklearn.decomposition import PCA


class herdingspikes(object):
    def __init__(self, classpath=None, datapath=None):
        self.datapath = datapath
        if classpath is not None:
            #  we want this to load a preexisting instantiation
            #  of detected stuff
            raise NotImplementedError()

    def LoadDetected(self):
        sp = np.loadtxt('ProcessedSpikes')
        self.spikes = pd.DataFrame({'ch': sp[:, 0].astype(int),
                                    't': sp[:, 1].astype(int),
                                    'Amplitude': sp[:, 2],
                                    'x': sp[:, 3],
                                    'y': sp[:, 4],
                                    'Shape': list(sp[:, 5:])
                                    })
        # save neighbour matrix and positions for the chip geometry
        f = open('neighbormatrix', 'r')
        self.neighs = []
        for l in f.readlines():
            self.neighs.append(np.array(l[:-2].split(',')).astype(int))
        f.close()
        f = open('positions', 'r')
        self.pos = []
        for l in f.readlines():
            self.pos.append(np.array(l[:-2].split(',')).astype(int))
        f.close()
        self.IsClustered = False

    def DetectFromRaw(self, num_channels, num_recording_channels, spike_delay,
                      spike_peak_duration, noise_duration, noise_amp,
                      max_neighbors, to_localize, cutout_length, fps, threshold,
                      maa=0, maxsl=12, minsl=3, ahpthr=0):
        detectData(self.datapath, num_channels, num_recording_channels,
                   spike_delay, spike_peak_duration, noise_duration, noise_amp,
                   max_neighbors, to_localize, cutout_length, fps, threshold,
                   maa, maxsl, minsl, ahpthr)
        # reload data into memory
        self.LoadDetected()

    def PlotTracesChannels(self, eventid, ax=None):
        if ax is None:
            ax = plt.gca()
        pos, neighs = self.pos, self.neighs
        data = np.fromfile(self.datapath,
                           dtype=np.int16).reshape((1800000, 385))

        event = self.spikes.loc[eventid]
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

    def CombinedMeanShift(self, bandwidth, alpha,
                          pca_ncomponents=2, pca_whiten=True,
                          bin_seeding=True, min_bin_freq=10, n_jobs=1):
        pca = PCA(n_components=pca_ncomponents, whiten=pca_whiten)
        cutouts_pca = pca.fit_transform(np.array(list(self.spikes.Shape)))
        fourvec = np.vstack(([self.spikes.x], [self.spikes.y],
                             alpha*cutouts_pca.T)).T
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=bin_seeding,
                       min_bin_freq=min_bin_freq, n_jobs=n_jobs)
        ms.fit(fourvec)
        self.spikes.cl = ms.labels_
        self.NClusters = len(np.unique(ms.labels_))
        assert np.max(self.spikes.cl)+1 == self.NClusters
        assert len(ms.cluster_centers_) == self.NClusters
        print("Number of estimated clusters:", self.NClusters)
        dic_cls = {'ctr_x': ms.cluster_centers_[:, 0],
                   'ctr_y': ms.cluster_centers_[:, 1],
                   'Color': 1.*np.random.permutation(
                        self.NClusters)/self.NClusters}
        # in_cl = []
        # for i in np.arange(ms.cluster_centers_.shape[1]):
        #     in_cl.append(list(np.where(self.spikes.cl == i)[0]))
        # dic_cls['Content'] = in_cl
        # dic_cls['AvgAmpl'] = [np.mean(self.spikes.Amplitude[cl])
        #                       for cl in in_cl]
        # dic_cls['Size'] = [len(cl) for cl in in_cl]
        self.clusters = pd.DataFrame(dic_cls)
        self.IsClustered = True

    def PlotShapes(self, units, nshapes=100, ncols=4):
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
