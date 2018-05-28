import logging
import os
import pickle
from collections import OrderedDict
from contextlib import redirect_stdout
from multiprocessing import Process, Queue, cpu_count

import numpy as np
import psutil
from .hs2 import HSClustering
from skopt import gp_minimize
from skopt.space import Integer, Real
from tqdm import tqdm


class OptimiseParameters(object):
    """ This class provides a simple way to optimise detection and sorting parameters of a spike sorting toolkit. The problem of finding best parameters can be reframed as finding a minimum of such black-box function f(x) that given a vector of (hyper)parameters x it performs detection/sorting and returns some scalar utility value reflecting the detection/sorting correctness. This minimum is searched for using Bayesian Optimisation with a Gaussian Process to approximate f(x) and the Expected Improvement acquisition function to select next values of x to try. Together with parameter regularisation, this provides a powerful framework for both optimisation and validation of spike sorting toolkits.

    To automatically evaluate the correctness of detection/sorting a paired ground-truth(GT) recording is required: a spiketrain with timestamps of true events.

    The correctness (=f(x)) of the detection performance is then determined as the number of detections missed on a particular channel (i.e. FNs). For sorting, we examine the cluster with most GT spikes, define precision & recall using that cluster and then define the utility as a negative F_1 score.

    Usage:
        1. Create an OptimiseParameters object by calling its constructor with a GT spiketrain, a HSDetection object, a HSClustering object and a Probe object. Use optimise_detection and optimise_clustering flags to specify parameters of which part of the pipeline shall be optimised.
        2. Call run().

    Parameters:
        closest_ch - a channel on the probe closest to the ground-truth recording device. I.E. a channel on which we expect the detection to be perfect.

        detec_params_to_opt - dictionary of parameters to be optimised together with respective ranges of their values defined as {"parameter_name": (lowest_value_to_be_tried, highest_value_to_be_tried), "parameter2_name": ...}.

        clust_params_to_opt - see detec_params_to_opt.

        true_positive_timewindow - (miliseconds). when optimising detection, how long to look before and after a GT spike for a spike on a channel on the probe.

        detec_run_schedule - a tuple where the *second* element is the number of random initiations before the optimisation begins, and the *first* element is the number of iterations the optimisation should be run for, including the random initiations.

        clust_run_schedule - see detec_run_schedule.

        clust_max_value - the current implementation (Mar2018) of HSClustering throws an error when parameters are set to ridiculous values. When this happens, what utility value should be returned.

        detec_outfile, clust_outfile - paths where optimisation results will be stored in a pickled format.
    """

    def __init__(self, gt_spiketrain, closest_ch, Probe, HSDetection,
                 detec_params_to_opt, HSClustering, clust_params_to_opt,
                 optimise_detection=True, optimise_clustering=True,
                 true_positive_timewindow=0.3, detec_run_schedule=[150, 100],
                 clust_run_schedule=[150, 100], clust_max_value=0,
                 detec_outfile='result_optim_params_detect',
                 clust_outfile='result_optim_params_clust'):

        logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

        self.gt_spiketrain = gt_spiketrain
        self.closest_ch = closest_ch

        self.optimise_detection = optimise_detection
        if optimise_detection:
            self.HSD = HSDetection
            # OrderedDict guarantees the correspondence between
            # the order of .keys() and .values()
            self.detec_params = OrderedDict(detec_params_to_opt)
            self.detec_outfile = detec_outfile
            self.detec_run_schedule = detec_run_schedule
            self.true_positive_timewindow = true_positive_timewindow / 1000 * Probe.fps  # in frames

        self.optimise_clustering = optimise_clustering
        if optimise_clustering:
            self.HSC = HSClustering
            self.clust_params = OrderedDict(clust_params_to_opt)
            self.clust_outfile = clust_outfile
            self.clust_run_schedule = clust_run_schedule
            self.clust_max_value = clust_max_value
            self.clust_utility_values = []
            self.neighbour_chs = Probe.neighbors[self.closest_ch].tolist()

        # Double-check threshold parameter is included for detection oprimisation
        if optimise_detection:
            assert 'threshold' in self.detec_params, "The detec_params_to_opt dictionary doesn't contain an entry for 'threshold' parameter. This parameter is required as it is used for detection regularisation."

    def run(self):
        """ Runs the optimisation of detection or sorting, depending of the optimise_detection and optimise_clustering flags.

        For each, saves results in detec_outfile or clust_outfile respectively. If both detection and sorting are being optimised, the most optimal detection parameters found are used for optimising sorting.
        The most optimal clustering configuration is saved into Clusters.hdf5

        Returns None, unless sorting is optimised, in which case returns a new HSClustering object with the optimal clustering performance.
        """
        if self.optimise_detection:
            results = self.optimise(
                self.detec_params, self.detection_wrapper, self.detec_run_schedule)

            # Perform detection with best parameters found
            best_parameters = {key: results.x[i] for i, key in enumerate(self.detec_params)}
            self.detection(best_parameters)
            TPs, FNs = self.detection_evaluate()

            # Add metadata and metrics to the results object
            results['TPs'] = TPs
            results['FNs'] = FNs
            results['n_GT_spikes'] = len(self.gt_spiketrain)
            results['detection_params'] = best_parameters

            # Hack to remove reference to a non-pickable function
            results.specs['args']['func'] = None

            self.save_results(results, self.detec_outfile)

        if self.optimise_clustering:
            # If both detection and clustering are being optimised, use best results
            # from detection optimisation
            if self.optimise_detection:
                self.HSC = HSClustering(self.HSD)
                # only GT spikes detected on the probe (TPs) will be used as GT for clustering
                self.clust_GT = self.load_results(self.detec_outfile)['TPs'][self.closest_ch]
            else:
                self.clust_GT = self.gt_spiketrain

            results = self.optimise(
                self.clust_params, self.clustering_wrapper, self.clust_run_schedule)

            # As the last step, cluster with best parameters and save the HSClustering object.
            _ = self.clustering_wrapper(results.x)
            self.HSC.SaveHDF5('Clusters.hdf5')

            # Add more metadata to the results object
            results['utility_values'] = self.clust_utility_values
            results['clustering_parameters'] = {key: results.x[i]
                                                for i, key in enumerate(self.clust_params)}
            results['TPs'], results['FPs'], results['most_popular_cluster'] = self.clustering_evaluate()
            results['Ps'] = self.clust_GT

            # Hack to remove reference to a non-pickable function
            results.specs['args']['func'] = None

            self.save_results(results, self.clust_outfile)

            return self.HSC

    def optimise(self, parameter_definitions, function, run_schedule):
        # Parse parameter definitions to a list of skopt dimensions
        dimensions = []
        for (name, rang) in parameter_definitions.items():
            if type(rang[1]) is int:
                dimensions.append(Integer(low=rang[0], high=rang[1], name=name))
            else:
                dimensions.append(Real(low=rang[0], high=rang[1], name=name))

        # Bayesian Optimisation using gaussian process(GP) using expected
        # improvement (EI) as an acquisition function
        results_object = gp_minimize(function,
                                     dimensions,
                                     acq_func='EI',
                                     noise=1e-10,
                                     n_calls=run_schedule[0],
                                     n_random_starts=run_schedule[1],
                                     n_jobs=-1,
                                     verbose=True
                                     )

        return results_object

    def save_results(self, obj, file_name):
        # Pickle the object into a new file
        with open("{}.pickle".format(file_name), 'wb+') as f:
            pickle.dump(obj, f)

    def load_results(self, file_name):
        with open("{}.pickle".format(file_name), 'rb') as f:
            return pickle.load(f)

    def detection_wrapper(self, chosen_values):
        """ A routine, which will run detection, evaluate the results and return a scalar value
        """
        # Label values with parameter names
        chosen_parameters = {}
        for i, key in enumerate(self.detec_params):
            chosen_parameters[key] = chosen_values[i]

        self.detection(chosen_parameters)

        TPs, FNs = self.detection_evaluate()
        missed_detection_per_channel = list(map(len, FNs))

        return self.detection_utility(missed_detection_per_channel, chosen_parameters['threshold'])

    def detection(self, parameters):
        # Hack to flush skopts verbose printing in jupyter notebook.
        print('', end='', flush=True)

        logging.info("Detecting spikes with parameters: {}".format(parameters))

        # Little hack to redirect the C++ output of HSDetection, parse it and log
        # it in a nice progress bar
        class ParseSTDOut:
            def __init__(self, pbar):
                self.pbar = pbar

            def write(self, line):
                if "# Analysing" in line and " frames;" in line:
                    no_frames = int(line.split()[2])
                    self.pbar.update(no_frames)

            def flush(self):
                pass

        with tqdm(desc='Detection', unit=' frames', unit_scale=1, total=self.HSD.probe.nFrames) as pbar:
            with redirect_stdout(ParseSTDOut(pbar)):
                self.HSD.SetAddParameters(parameters)
                self.HSD.DetectFromRaw(load=True)

    def detection_utility(self, missed_detection_per_channel, threshold_value):
        """ Defines the utility for evaluating detection performance.
        How many detections were missed on the most verbose channel REGULARISED by the threshold (we want to have threshold as strict as possible).
        """
        return missed_detection_per_channel[self.closest_ch] - threshold_value

    def detection_evaluate_per_channel(self, ch_spikes, neigh_spikes, channel, eval_queue):
        """ For every groundtruth spike, finds the closest spike from the current spiketrain and computes it's difference. If this difference is smaller than the true_positive_timewindow, we count it as a True Positive, else a True Negative. This procedure is first performed on the central channel spiketrain (ch_spikes) and then on the combined spiketrain from all neighbouring channels (neigh_spikes)
        """

        # If no spikes detected, count all gt_spikes as FNs
        if ch_spikes.shape[0] < 1 and neigh_spikes.shape[0] < 1:
            eval_queue.put((channel, [], self.gt_spiketrain))
            return

        # Helper method
        def update_bounds(idx):
            return (self.gt_spiketrain[idx] - self.true_positive_timewindow,
                    self.gt_spiketrain[idx] + self.true_positive_timewindow)

        idx = 0
        TPs, FNs, detected_gt_spikes = [], [], []
        lower_bound, upper_bound = update_bounds(idx)

        # First traverse the spiketrain of the central channel
        for raw_spike in ch_spikes:
            if raw_spike >= lower_bound:
                # If raw_spikes 'ahead' of GT spikes, advance gt_spikes to appropriate level
                while raw_spike > upper_bound:
                    idx += 1
                    if idx >= self.gt_spiketrain.shape[0]:
                        break
                    lower_bound, upper_bound = update_bounds(idx)

                if idx >= self.gt_spiketrain.shape[0]:
                    break

                # If raw_spike in +- true_positive_timewindow from gt_spike
                if lower_bound <= raw_spike and raw_spike <= upper_bound:
                    TPs.append(raw_spike)
                    detected_gt_spikes.append(idx)
                    idx += 1
                    if idx >= self.gt_spiketrain.shape[0]:
                        break
                    lower_bound, upper_bound = update_bounds(idx)

        # Move the counter to first gt_spike not yet detected
        idx = 0
        while idx in detected_gt_spikes:
            idx += 1

        # If all gt spikes were detected on the central channel, return
        if idx >= self.gt_spiketrain.shape[0]:
            eval_queue.put((channel, sorted(TPs), sorted(FNs)))
            return

        lower_bound, upper_bound = update_bounds(idx)

        # Traverse the neighbour channels' spiketrain
        for raw_spike in neigh_spikes:
            if raw_spike >= lower_bound:
                # If raw_spikes 'ahead' of GT spikes, advance gt_spikes to appropriate level
                while raw_spike > upper_bound:
                    if idx not in detected_gt_spikes:
                        FNs.append(raw_spike)

                    idx += 1
                    if idx >= self.gt_spiketrain.shape[0]:
                        break
                    lower_bound, upper_bound = update_bounds(idx)

                if idx >= self.gt_spiketrain.shape[0]:
                    break

                if lower_bound <= raw_spike and raw_spike <= upper_bound:
                    if idx not in detected_gt_spikes:
                        TPs.append(raw_spike)
                    idx += 1
                    if idx >= self.gt_spiketrain.shape[0]:
                        break
                    lower_bound, upper_bound = update_bounds(idx)

        # If there are no more raw_spikes to match, but still GT_spikes to be matched.
        while idx < self.gt_spiketrain.shape[0]:
            FNs.append(self.gt_spiketrain[idx])
            idx += 1
            while idx in detected_gt_spikes:
                idx += 1

        eval_queue.put((channel, sorted(TPs), sorted(FNs)))

    def detection_evaluate(self, n_CPUs=cpu_count() - 1):
        """ Parallelised detection_evaluate_per_channel() over n_CPUs processes. Each process is given spiketrain from one channel (data-parallelism).
        """
        num_channels = self.HSD.probe.num_channels

        # Shared queue used by the parallel evaluation
        eval_queue = Queue()

        # Initialize some datastructures (since we care about ordering)
        TPs = [[] for _ in range(num_channels)]
        FNs = [[] for _ in range(num_channels)]
        ch = 0
        with tqdm(desc='Evaluation', unit=' channels', total=num_channels) as pbar:
            while ch < num_channels:
                processes = []
                for _ in range(n_CPUs):
                    if ch == num_channels:
                        break

                    # Collate timestamps of all spikes on all neighbouring channels together
                    neighs_only = self.HSD.probe.neighbors[ch].tolist()
                    neighs_only.remove(ch)
                    neigh_spiketrain = self.HSD.spikes.loc[self.HSD.spikes.ch.isin(
                        neighs_only)]['t'].values
                    # Collate timestamps of all spikes on the central channel
                    ch_spiketrain = self.HSD.spikes.loc[self.HSD.spikes.ch == ch]['t'].values
                    # Make sure we don't run out of memory by spawning one more processes
                    if psutil.Process(os.getpid()).memory_info().rss > psutil.virtual_memory().available:
                        break

                    # Create a process, which evaluates TPs on all neighbouring channels of ch
                    processes.append(Process(target=self.detection_evaluate_per_channel,
                                             args=[ch_spiketrain, neigh_spiketrain, ch, eval_queue]))
                    ch += 1
                    processes[-1].start()

                # Collates results into appropriate lists. NOTE: Python is funny and this
                # also auto-joins child processes, as they are auto-joined as soon as
                # their contribution to Queue is consumed.
                for _ in range(len(processes)):
                    channel, tp, tn = eval_queue.get()
                    TPs[channel], FNs[channel] = tp, tn
                    pbar.update(1)

        eval_queue.close()

        return TPs, FNs

    def clustering_wrapper(self, chosen_values):
        clustering_params = {key: chosen_values[i] for i, key in enumerate(self.clust_params)}
        logging.info("Clustering spikes with parameters: {}".format(clustering_params))
        try:
            self.HSC.ShapePCA(
                pca_ncomponents=clustering_params['pca_ncomponents'], pca_whiten=True)
            self.HSC.CombinedClustering(alpha=clustering_params['alpha'],
                                        bandwidth=clustering_params['bandwidth'],
                                        bin_seeding=True,
                                        min_bin_freq=2,
                                        n_jobs=-1)
        # If the bandwidth is set too low, the algorithm throws an error.
        except ValueError:
            return self.clust_max_value

        TPs, FPs, most_popular_cluster = self.clustering_evaluate()
        utility = self.clustering_utility(TPs, FPs)

        self.clust_utility_values.append((utility, TPs.shape[0], FPs.shape[0]))

        return utility

    def clustering_utility(self, TPs, FPs):
        """ Defines the utility for optimising clustering performance.
        The F1 score, where in terms of GT and clustering:
        recall = how many GT spikes are in the most popular clusters
        precision = how many of the spikes in most popular cluster are GT spikes
        Most popular cluster is defined as the cluster containing most GT spikes.
        """
        if TPs.shape[0] == 0:
            return 0

        precision = TPs.shape[0] / (TPs.shape[0] + FPs.shape[0])
        recall = TPs.shape[0] / len(self.clust_GT)

        # Negated F1 score (since we are minimizing)
        return -1 * 2 * (precision * recall) / (precision + recall)

    def clustering_evaluate(self):
        """ Finds the cluster which contains the most GT spikes and identifies TPs and FPs:
         TP = a GT spike in this clusters
         FP = a spike in this cluster which is not in the GT
        """
        mask = self.HSC.spikes['t'].isin(
            self.clust_GT) & self.HSC.spikes['ch'].isin(self.neighbour_chs)
        Ps_in_neighbourhood = self.HSC.spikes.loc[mask]

        all_clusters = Ps_in_neighbourhood['cl'].values
        most_popular_cluster = np.bincount(all_clusters).argmax()

        TPs = Ps_in_neighbourhood.loc[all_clusters == most_popular_cluster]

        mask = ~self.HSC.spikes['t'].isin(self.clust_GT) & (
            self.HSC.spikes['cl'] == most_popular_cluster)
        FPs = self.HSC.spikes.loc[mask]

        return TPs, FPs, most_popular_cluster
