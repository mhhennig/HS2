from __future__ import division
from __future__ import absolute_import
import pandas as pd
import numpy as np
import h5py
import os
from pathlib import Path

from .detection_localisation.detect import detectData
from .detection_lightning.detect import HSDetectionLightning as detectDataLightning
from matplotlib import pyplot as plt
from .clustering.mean_shift_ import MeanShift
from sklearn.decomposition import PCA
from os.path import splitext
import warnings
from numpy.lib.stride_tricks import as_strided


def min_func(x):
    return x.min()


def max_func(x):
    return x.max()


class HSDetection(object):
    """
    This class provides a simple interface to the detection, localisation of
    spike data from dense multielectrode arrays according to the methods
    described in the following papers:

    Muthmann, J. O., Amin, H., Sernagor, E., Maccione, A., Panas, D.,
    Berdondini, L., ... & Hennig, M. H. (2015). Spike detection for large neural
    populations using high density multielectrode arrays. Frontiers in
    neuroinformatics, 9.

    Hilgen, G., Sorbaro, M., Pirmoradian, S., Muthmann, J. O., Kepiro, I. E.,
    Ullo, S., ... & Hennig, M. H. (2017). Unsupervised spike sorting for
    large-scale, high-density multielectrode arrays. Cell reports, 18(10),
    2521-2532.

    This is the original legacy implementation, here for backwards compatibility.
    This is now replaced by HSDetectionLightning.

    Usage:
        1. Create a HSDetection object by calling its constructor with a
    Probe object and all the detection parameters (see documentation there).
        2. Call DetectFromRaw.
        3. Save the result, or create a HSClustering object.

    Parameters
    ----------
    probe -- probe object, with a link to raw data
    to_localize : bool
        Set to False if spikes should only be detected, not localised (will break sorting)
    cutout_start : int
        Frame-based version of left_cutout_ms
    cutout_end : int
        Deprecated, frame-based version of right_cutout_ms
    threshold : int
        Detection threshold
    maa : int
        Minimum average amplitude
    maxsl : int
        Deprecated, frame-based version of spk_evaluation_time
    minsl : int
        Deprecated, frame-based version of amp_evaluation_time
    ahpthr : int
        Threshold for the repolarisation.
    out_file_name : str
        base file name (without extension) for the output files
    file_directory_name : str
        directory where output is saved
    save_all : bool
        Whether to stor all detection info (for debug)
    left_cutout_ms : float
        Milliseconds before the spike to be included in the spike shape
    right_cutout_ms : float
        Milliseconds after the spike to be included in the spike shape
    amp_evaluation_time : float
        Milliseconds during which the trace is integrated, for the purposed
        of evaluating amplitude, used for later comparison with 'maa'
    spk_evaluation_time : float
        Dead time in ms after spike peak, used for further triaging of spike shape
    """

    def __init__(
        self,
        probe,
        to_localize=True,
        num_com_centers=1,
        cutout_start=None,
        cutout_end=None,
        threshold=20,
        maa=0,
        maxsl=None,
        minsl=None,
        ahpthr=0,
        out_file_name="ProcessedSpikes",
        file_directory_name="HS2_results",
        decay_filtering=False,
        save_all=False,
        left_cutout_time=1.0,
        right_cutout_time=2.2,
        amp_evaluation_time=0.4,
        spk_evaluation_time=1.7,
    ):
        self.probe = probe

        self.cutout_start = self._deprecate_or_convert(
            cutout_start, left_cutout_time, "cutout_start", "left_cutout_time"
        )
        self.cutout_end = self._deprecate_or_convert(
            cutout_end, right_cutout_time, "cutout_end", "right_cutout_time"
        )
        self.minsl = self._deprecate_or_convert(
            minsl, amp_evaluation_time, "minsl", "amp_evaluation_time"
        )
        self.maxsl = self._deprecate_or_convert(
            maxsl, spk_evaluation_time, "maxsl", "spk_evaluation_time"
        )

        self.to_localize = to_localize
        self.cutout_length = self.cutout_start + self.cutout_end + 1
        self.threshold = threshold
        self.maa = maa
        self.ahpthr = ahpthr
        self.decay_filtering = decay_filtering
        self.num_com_centers = num_com_centers
        self.sp_flat = None
        self.spikes = None

        # Make directory for results if it doesn't exist
        os.makedirs(file_directory_name, exist_ok=True)

        if out_file_name[-4:] == ".bin":
            file_path = os.path.join(file_directory_name, out_file_name)
            self.out_file_name = file_path
        else:
            file_path = os.path.join(file_directory_name, out_file_name)
            self.out_file_name = file_path + ".bin"
        self.save_all = save_all

    def _deprecate_or_convert(self, old_var, new_var, old_name, new_name):
        if old_var is not None:
            msg = (
                "{} is deprecated and will be removed. Set {} instead "
                "(in milliseconds). {} takes priority over {}!"
            )
            warnings.warn(msg.format(old_name, new_name, old_name, new_name))
            return int(old_var)
        else:
            return int(new_var * self.probe.fps / 1000 + 0.5)

    def SetAddParameters(self, dict_of_new_parameters):
        """
        Adds and merges dict_of_new_parameters with the current fields of the
        object. Uses the PEP448 convention to group two dics together.

        Parameters
        ----------
        dict_of_new_parameters : dict
            The dictionary of parameters to be updated.
        """
        self.__dict__.update(dict_of_new_parameters)

    def LoadDetected(self, file_name=None):
        """
        Reads a binary file with spikes detected with the DetectFromRaw()
        method.

        Parameters
        ----------
        file_name : str
            The name of the .bin file. Defaults to self.out_file_name.
        """
        if file_name is None:
            file_name = self.out_file_name

        if os.stat(file_name).st_size == 0:
            shapecache = np.asarray([]).reshape(0, 5)
            warnings.warn(
                "Loading an empty file {} . This usually happens when no spikes were"
                "detected due to the detection parameters being set too "
                "strictly".format(file_name)
            )
        else:
            #             self.sp_flat = np.memmap(file_name, dtype=np.int16, mode="r")
            self.sp_flat = np.memmap(file_name, dtype=np.int32, mode="r")
            assert (
                self.sp_flat.shape[0] // self.cutout_length + 5
                is not self.sp_flat.shape[0] / self.cutout_length + 5
            ), "spike data has wrong dimensions"  # ???
            shapecache = self.sp_flat.reshape((-1, self.cutout_length + 5))

        self.spikes = pd.DataFrame(
            {
                "ch": shapecache[:, 0],
                "t": shapecache[:, 1],
                "Amplitude": shapecache[:, 2],
                "x": shapecache[:, 3] / 1000,
                "y": shapecache[:, 4] / 1000,
                "Shape": list(shapecache[:, 5:]),
            },
            copy=False,
        )
        self.IsClustered = False
        print("Loaded " + str(self.spikes.shape[0]) + " spikes.")

    def DetectFromRaw(
        self, load=True, nFrames=None, tInc=100000, recording_duration=None
    ):
        """
        This function is a wrapper of the C function `detectData`. It takes
        the raw data file, performs detection and localisation, saves the result
        to HSDetection.out_file_name and loads the latter into memory by calling
        LoadDetected if load=True.

        Parameters
        ----------
        load : bool
            Load the detected spikes into memory when finished.
        nFrames : int
            Frame-based version of recording_duration.
        tInc : ibt
            Size of chunks to be read.
        recording_duration : int
            maximum time limit of recording (the rest will be ignored)
        """

        if nFrames is not None:
            warnings.warn(
                "nFrames is deprecated and will be removed. Leave "
                "this out if you want to read the whole recording, "
                "or set max_duration to set the limit (in seconds)."
            )
        elif recording_duration is not None:
            nFrames = int(recording_duration * self.probe.fps)

        detectData(
            probe=self.probe,
            file_name=str.encode(self.out_file_name[:-4]),
            to_localize=self.to_localize,
            sf=self.probe.fps,
            thres=self.threshold,
            cutout_start=self.cutout_start,
            cutout_end=self.cutout_end,
            maa=self.maa,
            maxsl=self.maxsl,
            minsl=self.minsl,
            ahpthr=self.ahpthr,
            num_com_centers=self.num_com_centers,
            decay_filtering=self.decay_filtering,
            verbose=self.save_all,
            nFrames=nFrames,
            tInc=tInc,
        )

        if load:
            self.LoadDetected()

    def PlotTracesChannels(
        self,
        eventid,
        ax=None,
        window_size=100,
        show_channels=True,
        ascale=1.0,
        show_channel_numbers=True,
        show_loc=True,
    ):
        """
        Draw a figure with an electrode and its neighbours, showing the raw
        traces and events. Note that this requires loading the raw data in
        memory again.

        Parameters
        ----------
        eventid : int
            Centers, spatially and temporally, the plot to a specific event id.
        ax : matplotlib axes
            A matplotlib axes object where to draw. Defaults to current axis.
        window_size : int
            Number of samples shown around a spike.
        show_channels : bool
            Show bubbles corresponding to electrode locations.
        ascale : float
            Scaling factor for traces.
        show_channel_numbers : bool
            Whether to print electrode numbers next to them.
        show_loc : bool
            Whether to show a red point where the spike was localised.
        """
        if ax is None:
            ax = plt.gca()

        pos, neighs = self.probe.positions, self.probe.neighbors

        event = self.spikes.loc[eventid]
        print("Spike detected at channel: ", event.ch)
        print("Spike detected at frame: ", event.t)
        print("Spike localised in position", event.x, event.y)
        cutlen = len(event.Shape)

        # compute distance between electrodes, for scaling
        distances = np.abs(pos[event.ch][0] - pos[neighs[event.ch]][:, 0])
        interdistance = np.min(distances[distances > 0])
        scale = interdistance / 110.0 * ascale

        # scatter of the large grey balls for electrode location
        x = pos[(neighs[event.ch], 0)]
        y = pos[(neighs[event.ch], 1)]
        if show_channels:
            plt.scatter(x, y, s=1600, alpha=0.2)

        ws = window_size // 2
        ws = max(ws, 1 + self.cutout_start, 1 + self.cutout_end)
        t1 = np.max((0, event.t - ws))
        t2 = event.t + ws

        trange = (np.arange(t1, t2) - event.t) * scale
        start_bluered = event.t - t1 - self.cutout_start
        trange_bluered = trange[start_bluered : start_bluered + cutlen]

        data = self.probe.Read(t1, t2).reshape((t2 - t1, self.probe.num_channels))
        # remove offsets
        data = data - data[0]

        # grey and blue traces
        for i, n in enumerate(neighs[event.ch]):
            dist_from_max = np.sqrt(
                (pos[n][0] - pos[event.ch][0]) ** 2
                + (pos[n][1] - pos[event.ch][1]) ** 2
            )
            if n in self.probe.masked_channels:
                col = "g"
            elif dist_from_max <= self.probe.inner_radius:
                col = "orange"
            else:
                col = "b"
            plt.plot(pos[n][0] + trange, pos[n][1] + data[:, n] * scale, "gray")
            plt.plot(
                pos[n][0] + trange_bluered,
                pos[n][1] + data[start_bluered : start_bluered + cutlen, n] * scale,
                col,
            )

        # red overlay for central channel
        plt.plot(
            pos[event.ch][0] + trange_bluered,
            pos[event.ch][1] + event.Shape * scale,
            "r",
        )
        inner_radius_circle = plt.Circle(
            (pos[event.ch][0], pos[event.ch][1]),
            self.probe.inner_radius,
            color="red",
            fill=False,
        )
        ax.add_artist(inner_radius_circle)

        # red dot of event location
        if show_loc:
            plt.scatter(event.x, event.y, s=80, c="r")

        # electrode numbers
        if show_channel_numbers:
            for i, txt in enumerate(neighs[event.ch]):
                ax.annotate(txt, (x[i], y[i]))
        ax.set_aspect("equal")
        return ax

    def PlotDensity(self, binsize=1.0, invert=False, ax=None):
        """
        Plot the density of the spikes on the probe.

        Parameters
        ----------
        binsize : float
            The size of the bins for the density plot
        invert : bool
            Invert the x and y axes
        ax : matplotlib axis
            The axis to plot the density

        Returns
        -------
        h : array
            The histogram of the density
        xb : array
            The x bins
        yb : array
            The y bins
        """

        if ax is None:
            ax = plt.gca()
        x, y = self.spikes.x, self.spikes.y
        if invert:
            x, y = y, x
        binsx = np.arange(x.min(), x.max(), binsize)
        binsy = np.arange(y.min(), y.max(), binsize)
        h, xb, yb = np.histogram2d(x, y, bins=[binsx, binsy])
        ax.imshow(
            np.clip(np.log10(h), 1e-10, None),
            extent=[xb.min(), xb.max(), yb.min(), yb.max()],
            interpolation="none",
            origin="lower",
        )
        return h, xb, yb

    def PlotAll(self, invert=False, ax=None, max_show=100000, **kwargs):
        """
        Plot spikes currently stored in the class, in (x, y) space.

        Parameters
        ----------
        invert : bool
            Invert the x and y axes
        ax : matplotlib axis
            The axis to plot the spikes
        max_show : int
            The maximum number of spikes to show
        **kwargs : dict
            Additional arguments for the scatter plot

        Returns
        -------
        ax : matplotlib axis
            The axis with the spikes
        """
        if ax is None:
            ax = plt.gca()
        x, y = self.spikes.x, self.spikes.y
        if invert:
            x, y = y, x
        if self.spikes.shape[0] > max_show:
            inds = np.random.choice(self.spikes.shape[0], max_show, replace=False)
            print("We have", self.spikes.shape[0], "spikes, only showing", max_show)
        else:
            inds = np.arange(self.spikes.shape[0])
        ax.scatter(x[inds], y[inds], **kwargs)
        return ax

    def Cluster(self):
        """
        Create a HSClustering object.

        Returns
        -------
        HSClustering
            The HSClustering object
        """
        return HSClustering(self)


class HSDetectionLightning(object):
    """
    This class provides a simple interface to the detection, localisation of spikes
    from dense multielectrode arrays. The detection is based on the methods described
    in the following papers:

    Muthmann, J. O., Amin, H., Sernagor, E., Maccione, A., Panas, D.,
    Berdondini, L., ... & Hennig, M. H. (2015). Spike detection for large neural
    populations using high density multielectrode arrays. Frontiers in
    neuroinformatics, 9.

    Hilgen, G., Sorbaro, M., Pirmoradian, S., Muthmann, J. O., Kepiro, I. E.,
    Ullo, S., ... & Hennig, M. H. (2017). Unsupervised spike sorting for
    large-scale, high-density multielectrode arrays. Cell reports, 18(10),
    2521-2532.

    Usage:
        1. Create a HSDetection object by calling its constructor with a SpikeInterface
        Recording extractor and an optional detection parameter dictionary (see documentation).
        2. Call DetectFromRaw.
        3. Save the result, or create a HSClustering object.

    Parameters
    ----------
    rec : RecordingExtractor
        The recording extractor object
    params : dict
        The parameters for the spike detection. If None defaults will be used.
    """

    def __init__(self, rec, params=None):
        self.sp_flat = None
        self.spikes = None
        self.recording = rec
        self.num_segments = rec.get_num_segments()
        self.params = params
        if self.params == None:
            self.params = detectDataLightning.DEFAULT_PARAMS
        self.out_file_name = self.params["out_file"]
        if self.params["chunk_size"] == None:
            self.params["chunk_size"] = min((1e9 / rec.get_num_channels()), 1e6)
        if self.out_file_name is not None:
            out_dir = os.path.dirname(self.out_file_name)
            self.spikes_file = self.out_file_name + ".hdf5"

    def DetectFromRaw(self):
        """
        Detect spikes from the raw data.
        """
        det = detectDataLightning(self.recording, self.params)
        sp = det.detect()
        self.spikes = pd.DataFrame(
            {
                "ch": sp[0]["channel_ind"],
                "t": sp[0]["sample_ind"],
                "Amplitude": sp[0]["amplitude"],
                "x": sp[0]["location"][:, 0],
                "y": sp[0]["location"][:, 1],
                "Shape": list(-sp[0]["spike_shape"]),
            },
            copy=False,
        )
        # check if any spikes
        assert self.spikes.shape[0] > 0, "No spikes detected"
        # write spikes dict to hdf5 (without shapes)
        if self.out_file_name is not None:
            print(f"writing spikes to {self.spikes_file}")
            h = h5py.File(self.spikes_file, "w")
            for key in ["ch", "t", "Amplitude", "x", "y"]:
                h.create_dataset(key, data=self.spikes[key])
            h["cutout_length"] = self.spikes["Shape"][0].shape[0]
            h.close()

    def PlotTracesChannels(
        self,
        eventid,
        ax=None,
        window_size=100,
        show_channels=True,
        ascale=1.0,
        show_channel_numbers=True,
        show_loc=True,
    ):
        print("not implemented")

    def PlotDensity(self, binsize=1.0, invert=False, ax=None):
        """
        Plot the density of the spikes on the probe.

        Parameters
        ----------
        binsize : float
            The size of the bins for the density plot
        invert : bool
            Invert the x and y axes
        ax : matplotlib axis
            The axis to plot the density

        Returns
        -------
        h : array
            The histogram of the density
        xb : array
            The x bins
        yb : array
            The y bins
        """
        if ax is None:
            ax = plt.gca()
        x, y = self.spikes.x, self.spikes.y
        if invert:
            x, y = y, x
        binsx = np.arange(x.min(), x.max(), binsize)
        binsy = np.arange(y.min(), y.max(), binsize)
        h, xb, yb = np.histogram2d(x, y, bins=[binsx, binsy])
        ax.imshow(
            np.clip(np.log10(h), 1e-10, None),
            extent=[xb.min(), xb.max(), yb.min(), yb.max()],
            interpolation="none",
            origin="lower",
        )
        return h, xb, yb

    def PlotAll(self, invert=False, ax=None, max_show=100000, **kwargs):
        """
        Plot spikes currently stored in the class, in (x, y) space.

        Parameters
        ----------
        invert : bool
            Invert the x and y axes
        ax : matplotlib axis
            The axis to plot the spikes
        max_show : int
            The maximum number of spikes to show
        **kwargs : dict
            Additional arguments for the scatter plot

        Returns
        -------
        ax : matplotlib axis
            The axis with the spikes
        """
        if ax is None:
            ax = plt.gca()
        x, y = self.spikes.x, self.spikes.y
        if invert:
            x, y = y, x
        if self.spikes.shape[0] > max_show:
            inds = np.random.choice(self.spikes.shape[0], max_show, replace=False)
            print("We have", self.spikes.shape[0], "spikes, only showing", max_show)
        else:
            inds = np.arange(self.spikes.shape[0])
        ax.scatter(x[inds], y[inds], **kwargs)
        return ax

    def Cluster(self):
        """
        Create a HSClustering object.

        Returns
        -------
        HSClustering
            The HSClustering object
        """
        return HSClustering(self)


class HSClustering(object):
    """
    This class provides an easy interface to the clustering of spikes based
    on spike location on the chip and spike waveform, as described in:

    Hilgen, G., Sorbaro, M., Pirmoradian, S., Muthmann, J. O., Kepiro, I. E.,
    Ullo, S., ... & Hennig, M. H. (2017). Unsupervised spike sorting for
    large-scale, high-density multielectrode arrays. Cell reports, 18(10),
    2521-2532.

    Parameters
    ----------
    arg1 : str or list or HSDetection or HSDetectionLightning object
        Detected spikes. This can either be a HSDetection/HSDetectionLightning object,
        a file name for a list of file names. Lists will be concatenated into a single dataframe.
    cutout_length : int
        The length of the cutouts in the binary file, only required if a file or a
        list of files is passed as first argument and legacy==True.
    legacy : bool
        Whether to load the spikes from the legacy format (True) or the new Lightning format.
    **kwargs : dict
        Additional arguments to load the data.
    """

    def __init__(self, arg1, legacy=False, cutout_length=None, **kwargs):
        self.shapecache = []
        if type(arg1) == pd.core.frame.DataFrame:
            print("Reading spikes from Dataframe")
            self.spikes = arg1
            self.expinds = [0]
        elif type(arg1) == HSDetectionLightning or type(arg1) == HSDetection:
            print("Reading spikes from detection")
            self.spikes = arg1.spikes
            self.expinds = [0]
        else:
            if type(arg1) == str:
                arg1 = [arg1]
            if legacy:
                for i, f in enumerate(arg1):
                    filetype = splitext(f)[-1]
                    not_first_file = i > 0
                    if filetype == ".hdf5":
                        _f = h5py.File(f, "r")
                        if "shapes" in list(_f.keys()):
                            _f.close()
                            self.LoadHDF5(f, append=not_first_file, **kwargs)
                        elif "Shapes" in list(_f.keys()):
                            _f.close()
                            self.LoadHDF5_legacy_detected(
                                f, append=not_first_file, **kwargs
                            )
                    elif filetype == ".bin":
                        if cutout_length is None:
                            raise ValueError(
                                "You must pass cutout_length for .bin files."
                            )
                        self.LoadBin(f, cutout_length, append=not_first_file)
                    else:
                        raise IOError("File format unknown. Expected .hdf5 or .bin")
            else:
                for i, f in enumerate(arg1):
                    shape_file = Path(f).with_suffix(".bin")
                    shape_file = shape_file.with_stem(f"{shape_file.stem}-0")
                    spikes_file = Path(f).with_suffix(".hdf5")
                    print(
                        f"loading spikes from {spikes_file}, shapes from {shape_file}"
                    )
                    spikes_data = h5py.File(spikes_file, "r")
                    cutout_length = spikes_data["cutout_length"][()]
                    shapes = np.memmap(
                        str(shape_file), dtype=np.int16, mode="r"
                    ).reshape(-1, cutout_length)
                    self.shapecache.append(shapes)
                    spikes = pd.DataFrame(
                        {
                            "Shape": list(-self.shapecache[-1]),
                            "ch": spikes_data["ch"],
                            "t": spikes_data["t"],
                            "Amplitude": spikes_data["Amplitude"],
                            "x": spikes_data["x"],
                            "y": spikes_data["y"],
                        },
                        copy=False,
                    )
                    if i > 0:
                        self.expinds.append(len(self.spikes))
                        self.spikes = pd.concat(
                            [self.spikes, spikes], ignore_index=True
                        )
                        self.filelist.append(arg1)
                    else:
                        self.spikes = spikes
                        self.expinds = [0]
                        self.filelist = [arg1]

    def CombinedClustering(
        self, alpha, clustering_algorithm=None, cluster_subset=None, **kwargs
    ):
        """
        Clusters spikes based on their location and additional features (e.g.
        PCA projections of the spike shapes at their peak channel).
        Cluster memberships are available in the dataframe HSClustering.spikes (column
        cl), and cluster information in HSClustering.clusters.

        Parameters
        ----------
        alpha : float
            The weight given to the other features, relative to spatial components
            (which have weight 1.)
        clustering_algorithm : sklearn.cluster object
            The clustering algorithm to use. If None, a MeanShift object will be
            instantiated.
        cluster_subset : int
            Number of spikes used to build clusters, spikes are then assigned to the
            nearest by Euclidean distance
        **kwargs : dict
            Additional arguments are passed to the clustering class. This may include
            n_jobs > 1 for parallelisation.
        """
        try:
            fourvec = np.vstack(
                (self.spikes.x, self.spikes.y, alpha * self.features.T)
            ).T
        except AttributeError:
            fourvec = np.vstack((self.spikes.x, self.spikes.y)).T
            print("Warning: no PCA or other features available, location only!")

        print("Clustering...")
        if clustering_algorithm is None:
            clusterer = MeanShift(**kwargs)
        else:
            clusterer = clustering_algorithm(**kwargs)

        if cluster_subset is not None:
            print("Using", cluster_subset, "out of", self.spikes.shape[0], "spikes...")
            inds = np.sort(
                np.random.choice(
                    self.spikes.shape[0], int(cluster_subset), replace=False
                )
            )
            clusterer.fit(fourvec[inds])
            self.NClusters = len(np.unique(clusterer.labels_))
            print("Number of estimated units:", self.NClusters)
            print("Predicting cluster labels for", self.spikes.shape[0], "spikes...")
            self.spikes["cl"] = clusterer.predict(fourvec)
        else:
            print("Clustering " + str(self.spikes.shape[0]) + " spikes...")
            self.spikes["cl"] = clusterer.fit_predict(fourvec)
            self.NClusters = len(np.unique(self.spikes["cl"]))
            print("Number of estimated units:", self.NClusters)

        # methods like DBSCAN assign '-1' to unclustered data
        # here we replace these by a new cluster at the end of the list
        if self.spikes.cl.min() == -1:
            print(
                "There are",
                (self.spikes.cl == -1).sum(),
                "unclustered events, these are now in cluster number ",
                self.NClusters - 1,
            )
            self.spikes.loc[self.spikes.cl == -1, "cl"] = self.NClusters - 1

        _cl = self.spikes.groupby(["cl"])
        _x_mean = _cl.x.mean()
        _y_mean = _cl.y.mean()
        # this computes average amplitudes, disabled for now
        # _avgAmpl = _cl.min_amp.mean()
        _cls = _cl.cl.count()
        _color = 1.0 * np.random.permutation(self.NClusters) / self.NClusters
        dic_cls = {"ctr_x": _x_mean, "ctr_y": _y_mean, "Color": _color, "Size": _cls}
        #    'AvgAmpl': _avgAmpl
        # }
        self.clusters = pd.DataFrame(dic_cls)
        self.IsClustered = True

    def ShapePCA(
        self,
        pca_ncomponents=2,
        pca_whiten=True,
        chunk_size=1000000,
        custom_decomposition=None,
    ):
        """
        Computes the principal components of the spike shapes contained in the class,
        and computes projections of the spikes on these components. The features are
        saved to HSClustering.features, to be used for clustering.

        Parameters
        ----------
        pca_ncomponents : int
            Number of principal components to be used
        pca_whiten : bool
            Whiten data before PCA
        chunk_size : int
            Maximum number of shapes to be used to find PCs
        custom_decomposition : sklearn.decomposition object
            A custom instance of a sklearn decomposition object (such as instances PCA
            or FastICA), to be used for custom dimensionality reduction. pca_ncomponents
            and pca_whiten are ignored if provided.
        """

        def custom_roll(arr, r_tup):
            m = np.asarray(r_tup)
            arr_roll = arr[
                :, [*range(arr.shape[1]), *range(arr.shape[1] - 1)]
            ].copy()  # need `copy`
            strd_0, strd_1 = arr_roll.strides
            n = arr.shape[1]
            result = as_strided(arr_roll, (*arr.shape, n), (strd_0, strd_1, strd_1))
            return result[np.arange(arr.shape[0]), (n - m) % n]

        n_spikes = self.spikes.shape[0]
        if custom_decomposition is None:  # default is PCA
            _pca = PCA(n_components=pca_ncomponents, whiten=pca_whiten)
        else:  # Accepts an arbitrary sklearn.decomposition object instead of PCA
            _pca = custom_decomposition

        if n_spikes > chunk_size:
            print(
                f"Fitting dimensionality reduction using {chunk_size} out of {n_spikes} spikes..."
            )
            inds = np.sort(np.random.choice(n_spikes, chunk_size, replace=False))
            s = self.spikes.Shape.loc[inds].values.tolist()
        else:
            print("Fitting dimensionality reduction using all spikes...")
            s = self.spikes.Shape.values.tolist()
        # align peaks
        peak = np.argmin(np.diff(np.asarray(s).T, axis=0), axis=0)
        shift = np.min(peak) - peak
        s_roll = custom_roll(np.asarray(s), shift)

        _pca.fit(np.asarray(s_roll))

        print("...projecting...")
        _pcs = np.empty((n_spikes, pca_ncomponents))
        for i in range(n_spikes // chunk_size + 1):
            # is this the best way? Warning: Pandas slicing with .loc is different!
            s = self.spikes.Shape.loc[
                i * chunk_size : (i + 1) * chunk_size - 1
            ].values.tolist()
            # align peaks
            peak = np.argmin(np.diff(np.asarray(s).T, axis=0), axis=0)
            shift = np.min(peak) - peak
            s_roll = custom_roll(np.asarray(s), shift)
            _pcs[i * chunk_size : (i + 1) * chunk_size, :] = _pca.transform(s)

        self.pca = _pca
        self.features = _pcs
        print("...done")

    def _savesinglehdf5(
        self,
        filename,
        limits,
        compression,
        sampling,
        transpose=False,
        save_shapes=False,
    ):
        if limits is not None:
            spikes = self.spikes[limits[0] : limits[1]]
        else:
            spikes = self.spikes

        g = h5py.File(filename, "w")
        if transpose:
            g.create_dataset("data", data=np.vstack((spikes.y, spikes.x)))
        else:
            g.create_dataset("data", data=np.vstack((spikes.x, spikes.y)))
        if sampling is not None:
            g.create_dataset("Sampling", data=sampling)
        g.create_dataset("times", data=spikes.t)
        g.create_dataset("ch", data=spikes.ch)
        if self.IsClustered:
            if transpose:
                g.create_dataset("centres", data=self.clusters[["ctr_y", "ctr_x"]])
            else:
                g.create_dataset("centres", data=self.clusters[["ctr_x", "ctr_y"]])
            # g.create_dataset("centres", data=self.centers.T)
            g.create_dataset("cluster_id", data=spikes.cl)
        else:
            g.create_dataset("centres", data=[])
            g.create_dataset("cluster_id", data=[])

        g.create_dataset("exp_inds", data=self.expinds)
        # this is still a little slow (and perhaps memory intensive)
        # but I have not yet found a better way:
        # no longer save the spike shapes by default
        # if save_shapes and not spikes.empty:
        #     cutout_length = spikes.Shape.iloc[0].size
        #     sh_tmp = np.empty((cutout_length, spikes.Shape.size), dtype=int)
        #     for i, s in enumerate(spikes.Shape):
        #         sh_tmp[:, i] = s
        #     g.create_dataset("shapes", data=sh_tmp, compression=compression)
        # else:
        # g.create_dataset("shapes", data=[], compression=compression)
        g.close()

    def SaveHDF5(
        self,
        filename,
        compression=None,
        sampling=None,
        transpose=False,
        save_shapes=False,
    ):
        """
        Saves data, cluster centres and ClusterIDs to a hdf5 file. Offers
        compression of the shapes, 'lzf' appears a good trade-off between speed
        and performance.

        If filename is a single name, then all will be saved to a single file.
        If filename is a list of names of the same length as the number of
        experiments, one file per experiment will be saved.

        Parameters
        ----------
        filename : str or list
            The names of the file or list of files to be saved.
        compression : str
            Passed to HDF5, to save shapes compressed.
        sampling : float
            Provide sampling rate to include it in the file.
        transpose : bool
            Whether to swap x and y.  Default is False.
        save_shapes : bool
            Whether to save the spike shapes. Default is False.
        """
        if sampling is None:
            print(
                "# Warning: no sampling rate given, will be set to 0 in the hdf5 file."
            )
            sampling = 0

        if type(filename) == str:
            self._savesinglehdf5(filename, None, compression, sampling, transpose)
        elif type(filename) == list:
            if len(filename) != len(self.expinds):
                raise ValueError(
                    "Names list length does not correspond to "
                    "number of experiments in memory."
                )
            expinds = self.expinds + [len(self.spikes)]
            for i, f in enumerate(filename):
                self._savesinglehdf5(
                    f, [expinds[i], expinds[i + 1]], compression, sampling, transpose
                )
        else:
            raise ValueError("filename not understood")

    def LoadHDF5(self, filename, append=False, chunk_size=1000000, scale=1.0):
        """
        Load data, cluster centres and ClusterIDs from a hdf5 file created with
        legacy HS1 code folowing clustering.

        Parameters
        ----------
        filename : str
            The name of the file to load from.
        append : bool
            Append to data already in memory.
        chunk_size : int
            Read shapes in chunks of this size to avoid memory problems.
        scale : float
            Re-scale shapes by this factor (may be required for HS1 data).
        """
        g = h5py.File(filename, "r")
        print("Reading from clustered (HS1 or HS2) file " + filename)

        print(
            "Creating memmapped cache for shapes, reading in chunks of size",
            chunk_size,
            "and converting to integer...",
        )
        i = len(self.shapecache)
        self.shapecache.append(
            np.memmap(
                "tmp" + str(i) + ".bin",
                dtype=np.int32,
                mode="w+",
                shape=g["shapes"].shape[::-1],
            )
        )
        for i in range(g["shapes"].shape[1] // chunk_size + 1):
            tmp = (
                scale
                * np.transpose(g["shapes"][:, i * chunk_size : (i + 1) * chunk_size])
            ).astype(np.int32)
            inds = np.where(tmp > 20000)[0]
            tmp[inds] = 0
            print("Read chunk " + str(i + 1))
            if len(inds) > 0:
                print("Found", len(inds), "data points out of linear regime")
            print(len(self.shapecache), tmp.shape, i * chunk_size, (i + 1) * chunk_size)
            self.shapecache[-1][i * chunk_size : (i + 1) * chunk_size] = tmp[:]

        self.cutout_length = self.shapecache[-1].shape[1]
        print("Events: ", self.shapecache[-1].shape[0])
        print("Cut-out size: ", self.cutout_length)

        spikes = pd.DataFrame(
            {
                "ch": np.zeros(g["times"].shape[0], dtype=int),
                "t": g["times"],
                "Amplitude": np.zeros(g["times"].shape[0], dtype=int),
                "x": g["data"][0, :],
                "y": g["data"][1, :],
                "Shape": list(self.shapecache[-1]),
            },
            copy=False,
        )

        if "ch" in list(g.keys()):
            spikes.ch = g["ch"].value.T

        print("Getting spike amplitudes")
        spikes["min_amp"] = spikes.Shape.apply(min_func)
        spikes["Amplitude"] = spikes["min_amp"]

        if "centres" in list(g.keys()):
            self.centerz = g["centres"].value
            if len(self.centerz) < 5:
                print("WARNING Hack: Assuming HS1 data format")
                self.centerz = self.centerz.T
            self.NClusters = len(self.centerz)
            print("Number of clusters: ", self.NClusters)
            spikes["cl"] = g["cluster_id"]

            inds = spikes.groupby(["cl"]).cl.count().index
            _cl = spikes.groupby(["cl"])
            # this computes average amplitudes, disabled for now
            # _avgAmpl = np.zeros(self.NClusters)
            # _avgAmpl[inds] = _cl.min_amp.mean()
            _cls = np.zeros(self.NClusters)
            _cls[inds] = _cl.cl.count()

            dic_cls = {
                "ctr_x": self.centerz[:, 0],
                "ctr_y": self.centerz[:, 1],
                "Color": 1.0 * np.random.permutation(self.NClusters) / self.NClusters,
                "Size": _cls,
            }
            self.clusters = pd.DataFrame(dic_cls)
            self.IsClustered = True
        else:
            self.IsClustered = False

        g.close()

        if append:
            self.expinds.append(len(self.spikes))
            self.spikes = pd.concat([self.spikes, spikes], ignore_index=True)
            self.filelist.append(filename)
        else:
            self.spikes = spikes
            self.expinds = [0]
            self.filelist = [filename]

    def LoadHDF5_legacy_detected(
        self, filename, append=False, chunk_size=1000000, scale=1.0
    ):
        """
        Load data, cluster centres and ClusterIDs from a hdf5 file created with
        the HS1 detector.

        Parameters
        ----------
        filename : str
            The name of the file to load from.
        append : bool
            Append to data already in memory.
        chunk_size : int
            Read shapes in chunks of this size to avoid memory problems.
        scale : float
            Re-scale shapes by this factor (may be required for HS1 data).
        """
        g = h5py.File(filename, "r")
        print("Reading from unclustered HS1 file " + filename)
        if scale == 1:
            scale = -1.0 * g["Ascale"].value
        print(
            "Creating memmapped cache for shapes, reading in chunks of size",
            chunk_size,
            "and converting to integer...",
        )
        i = len(self.shapecache)
        self.shapecache.append(
            np.memmap(
                "tmp" + str(i) + ".bin",
                dtype=np.int32,
                mode="w+",
                shape=g["Shapes"].shape,
            )
        )
        for i in range(g["Shapes"].shape[0] // chunk_size + 1):
            tmp = (
                (scale * g["Shapes"][i * chunk_size : (i + 1) * chunk_size, :].T)
                .astype(np.int32)
                .T
            )
            inds = np.where(tmp > 20000)[0]
            tmp[inds] = 0
            print("Read chunk " + str(i + 1))
            if len(inds) > 0:
                print("Found", len(inds), "data points out of linear regime")
            self.shapecache[-1][i * chunk_size : (i + 1) * chunk_size] = tmp[:]

        self.cutout_length = self.shapecache[-1].shape[1]
        print("Events: ", self.shapecache[-1].shape[0])
        print("Cut-out size: ", self.cutout_length)

        spikes = pd.DataFrame(
            {
                "ch": np.zeros(g["Times"].shape[0], dtype=int),
                "t": g["Times"],
                "Amplitude": (-g["Amplitudes"][:] * scale).astype(int),
                "x": g["Locations"][:, 0],
                "y": g["Locations"][:, 1],
                "Shape": list(self.shapecache[-1]),
            },
            copy=False,
        )

        self.IsClustered = False

        g.close()

        if append:
            self.expinds.append(len(self.spikes))
            self.spikes = pd.concat([self.spikes, spikes], ignore_index=True)
            self.filelist.append(filename)
        else:
            self.spikes = spikes
            self.expinds = [0]
            self.filelist = [filename]

    def LoadBin(self, filename, cutout_length, append=False):
        """
        Reads a binary file with spikes detected with the legacy DetectFromRaw() method.

        Parameters
        ----------
        filename : str
            The name of the file
        cutout_length : int
            The length of each spike shape, in frames (this is necessary
            because the data is stored as a 1d array)
        append : bool
            Append to data already im memory
        """
        # 5 here are the non-shape data columns
        print("# loading", filename)
        self.shapecache.append(
            np.memmap(filename, dtype=np.int32, mode="r").reshape(
                (-1, cutout_length + 5)
            )
        )
        assert self.shapecache[-1].shape[0] // (
            cutout_length + 5
        ) is not self.shapecache[-1].shape[0] / (
            cutout_length + 5
        ), "spike data has wrong dimensions"  # ???
        spikes = pd.DataFrame(
            {
                "ch": self.shapecache[-1][:, 0],
                "t": self.shapecache[-1][:, 1],
                "Amplitude": self.shapecache[-1][:, 2],
                "x": self.shapecache[-1][:, 3] / 1000,
                "y": self.shapecache[-1][:, 4] / 1000,
                "Shape": list(self.shapecache[-1][:, 5:]),
            },
            copy=False,
        )
        self.IsClustered = False

        # this computes average amplitudes, disabled for now
        # spikes['min_amp'] = spikes.Shape.apply(min_func)

        if append:
            self.expinds.append(len(self.spikes))
            self.spikes = pd.concat([self.spikes, spikes], ignore_index=True)
            self.filelist.append(filename)
        else:
            self.spikes = spikes
            self.expinds = [0]
            self.filelist = [filename]

    def PlotShapes(
        self,
        units,
        ncols=4,
        ylim=None,
        max_shapes=100,
    ):
        """
        Plot a sample of the spike shapes contained in a given set of clusters
        and their average.

        Parameters
        ----------
        units : list
            A list of the cluster IDs to be considered.
        ncols : int
            The number of columns under which to distribute the plots.
        ylim : list
            Limits of the vertical axis of the plots. If None, try to figure
            them out.
        """
        nrows = int(np.ceil(len(units) / ncols))
        cutouts = self.spikes.Shape

        # all this is to determine suitable ylims TODO probe should provide
        yoff = 0
        if ylim is None:
            meanshape = np.mean(np.vstack(cutouts.loc[:1000].values), axis=0)
            yoff = -meanshape[0]
            maxy, miny = meanshape.max() + yoff, meanshape.min() + yoff
            varshape = np.var(cutouts.loc[:1000].values, axis=0)
            varmin = varshape[np.argmin(meanshape)]
            varmax = varshape[np.argmax(meanshape)]
            maxy += 4.0 * np.sqrt(varmax)
            miny -= 2.0 * np.sqrt(varmin)
            ylim = [miny, maxy]

        plt.figure(figsize=(3 * ncols, 3 * nrows))
        for i, cl in enumerate(units):
            inds = np.where(self.spikes.cl == cl)[0][:max_shapes]
            meanshape = np.mean(np.vstack(cutouts.loc[inds].values), axis=0)
            yoff = -meanshape[0]

            plt.subplot(nrows, ncols, i + 1)
            plt.plot(
                np.vstack(cutouts.loc[inds[:50]].values).T,
                color=(0.8, 0.8, 0.8),
                lw=0.8,
            )

            plt.plot(meanshape + yoff, c=plt.cm.hsv(self.clusters.Color[cl]), lw=4)
            plt.ylim(ylim)
            plt.title("Cluster " + str(cl))
        plt.tight_layout()

    def PlotAll(
        self,
        invert=False,
        show_labels=False,
        ax=None,
        max_show=200000,
        fontsize=16,
        **kwargs,
    ):
        """
        Plots all the spikes currently stored in the class, in (x, y) space.
        If clustering has been performed, each spike is coloured according to
        the cluster it belongs to.

        Parameters
        ----------
        invert : bool
            Invert the x and y axes
        show_labels : bool
            Annotate each cluster centre with its cluster ID
        ax : matplotlib axis
            The axis to plot the spikes
        max_show : int
            The maximum number of spikes to show
        fontsize : int
            The size of the font for the labels
        **kwargs : dict
            Additional arguments for the scatter plot

        Returns
        -------
        ax : matplotlib axis
            The axis with the plot
        """
        n_spikes = self.spikes.shape[0]
        if ax is None:
            ax = plt.gca()
        x, y = self.spikes.x, self.spikes.y
        if invert:
            x, y = y, x
        if self.spikes.shape[0] > max_show:
            inds = np.random.choice(n_spikes, max_show, replace=False)
            print("We have", n_spikes, "spikes, only showing ", max_show)
        else:
            inds = np.arange(n_spikes)

        c = plt.cm.hsv(self.clusters.Color[self.spikes.cl])
        ax.scatter(x[inds], y[inds], c=c[inds], **kwargs)
        if show_labels and self.IsClustered:
            ctr_x, ctr_y = self.clusters.ctr_x, self.clusters.ctr_y
            if invert:
                ctr_x, ctr_y = ctr_y, ctr_x
            for cl in range(self.NClusters):  # TODO why is this here
                if ~np.isnan(ctr_y[cl]):  # hack, why NaN positions in DBScan?
                    ax.annotate(str(cl), [ctr_x[cl], ctr_y[cl]], fontsize=fontsize)
                    # seems this is a problem when zooming with x/ylim
        return ax

    def PlotNeighbourhood(
        self,
        cl,
        radius=1.0,
        show_cluster_numbers=True,
        alpha=0.4,
        show_unclustered=False,
        max_shapes=1000,
        figsize=(8, 6),
    ):
        """
        Plot all units and spikes in the neighbourhood of cluster cl.

        Parameters
        ----------
        cl : int
            The ID of the cluster to be shown
        radius : float
            Spikes are shown for units this far away from the cluster centre
        show_cluster_numbers : bool
            Whether to print cluster labels
        alpha : float
            Transparency of spike points
        show_unclustered : bool
            Whether to show unclustered spikes (left by certain clustering algorithms)
        max_shapes : int
            Maximum number of shapes to be plotted
        figsize : tuple
            The size of the figure

        Returns
        -------
        ax : matplotlib axis
            The axis with the plot.
        """

        plt.figure(figsize=figsize)

        cx, cy = self.clusters["ctr_x"][cl], self.clusters["ctr_y"][cl]
        dists = np.sqrt(
            (cx - self.clusters["ctr_x"]) ** 2 + (cy - self.clusters["ctr_y"]) ** 2
        )
        clInds = np.where(dists < radius)[0]

        ax = []
        ax.append(
            plt.subplot2grid(
                (len(clInds) + 1, 4),
                (0, 0),
                rowspan=len(clInds) + 1,
                colspan=3,
                facecolor="k",
            )
        )
        for i in range(len(clInds) + 1):
            ax.append(plt.subplot2grid((len(clInds) + 1, 4), (i, 3), colspan=1))
            ax[i + 1].axis("off")
            if i > 0:
                ax[i].sharey(ax[i + 1])

        for i_cl, cl_t in enumerate(clInds):
            cx, cy = self.clusters["ctr_x"][cl_t], self.clusters["ctr_y"][cl_t]
            inds = np.where(self.spikes.cl == cl_t)[0][:max_shapes]
            x, y = self.spikes.x[inds], self.spikes.y[inds]
            ax[0].scatter(
                x, y, color=plt.cm.hsv(self.clusters["Color"][cl_t]), s=3, alpha=alpha
            )
            if show_cluster_numbers:
                ax[0].text(cx - 0.1, cy, str(cl_t), fontsize=16, color="w")
            ax[i_cl + 2].plot(
                np.vstack(self.spikes.loc[inds[:50]].Shape.values).T,
                color=(0.8, 0.8, 0.8),
                lw=0.8,
            )
            if len(inds) > 1:
                ax[i_cl + 2].plot(
                    np.mean(np.vstack(self.spikes.loc[inds].Shape.values), axis=0),
                    color=plt.cm.hsv(self.clusters["Color"][cl_t]),
                    lw=2,
                )

        ax[0].axis("equal")

        # show unclustered spikes (if any)
        if show_unclustered:
            cx, cy = self.clusters["ctr_x"][cl], self.clusters["ctr_y"][cl]
            inds = np.where(self.spikes.cl == self.NClusters - 1)[0][:max_shapes]
            x, y = self.spikes.x[inds].values, self.spikes.y[inds].values
            dists = np.sqrt((cx - x) ** 2 + (cy - y) ** 2)
            spInds = np.where(dists < radius)[0]
            if len(spInds):
                ax[0].scatter(x[spInds], y[spInds], c="w", s=3)
                for i in spInds[:20]:
                    ax[1].plot(self.spikes.Shape[i], color=(0.4, 0.4, 0.4))
                ax[1].plot(np.mean(self.spikes.Shape[spInds].values, axis=0), color="k")
        return ax
