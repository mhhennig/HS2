# distutils: language=c++
# cython: language_level=3
# cython: annotation_typing=False

import warnings
from pathlib import Path
from typing import Optional

import cython
import numpy as np
from numpy.typing import NDArray

from ..recording import RealArray, Recording
from . import DEFAULT_PARAMS as _DEFAULT_PARAMS
from . import Params

bool_t = cython.typedef(_bool)  # type: ignore
int32_t = cython.typedef(_int32_t)  # type: ignore
p_i32 = cython.typedef(cython.pointer(int32_t))  # type: ignore
single = cython.typedef(cython.float)  # type: ignore
p_single = cython.typedef(cython.p_float)  # type: ignore
vector_i32 = cython.typedef(vector[int32_t])  # type: ignore

RADIUS_EPS: float = cython.declare(single, 1e-3)  # type: ignore  # 1nm


@cython.cclass
class HSDetectionLightning(object):
    """The spike detection algorithm in Herding Spikes. This class provides an \
    interface connecting SpikeInterface and the C++ implementation of the \
    algorithm.

    Reference:
    1. J.-O. Muthmann, H. Amin, E. Sernagor, A. Maccione, D. Panas, L. Berdondini, \
    U. S. Bhalla, M. H. Hennig. "Spike Detection for Large Neural Populations Using \
    High Density Multielectrode Arrays". In: *Frontiers in Neuroinformatics* 9 (2015).
    2. https://github.com/mhhennig/HS2

    Args:
        `recording` (`BaseRecording`): A recording from SpikeInterface.
        `params` (`dict`): A dictionary of algorithmic parameters. All keys \
            must be present (use `HSDetection.DEFAULT_PARAMS` as a start), but \
            additional keys are accepted and ignored.

    Returns (of `HSDetection.detect()`):
        `list[dict[str, np.ndarray]]`: A list of dictionary containing the \
                detection results for each recording segment.
        - sample_ind: (n,) array of sample/frame indices of detected spikes
        - channel_ind: (n,) array of channel indices
        - amplitude: (n,) array of amplitudes
        - location (optional): (n,2) array of spike locations if localization on
        - spike_shape (optional): (n,l) array of spike shape if shape saved
    """

    DEFAULT_PARAMS: Params = _DEFAULT_PARAMS

    recording: Recording = cython.declare(object)  # type: ignore
    num_segments: int = cython.declare(int32_t)  # type: ignore
    num_frames: list[int] = cython.declare(vector_i32)  # type: ignore

    num_channels: int = cython.declare(int32_t)  # type: ignore
    chunk_size: int = cython.declare(int32_t)  # type: ignore
    left_margin: int = cython.declare(int32_t)  # type: ignore

    rescale: bool = cython.declare(bool_t)  # type: ignore
    scale: NDArray[np.single] = cython.declare(np.ndarray)  # type: ignore
    offset: NDArray[np.single] = cython.declare(np.ndarray)  # type: ignore

    median_reference: bool = cython.declare(bool_t)  # type: ignore
    average_reference: bool = cython.declare(bool_t)  # type: ignore

    spike_duration: int = cython.declare(int32_t)  # type: ignore
    amp_avg_duration: int = cython.declare(int32_t)  # type: ignore
    threshold: float = cython.declare(single)  # type: ignore
    min_avg_amp: float = cython.declare(single)  # type: ignore
    max_AHP_amp: float = cython.declare(single)  # type: ignore

    positions: NDArray[np.single] = cython.declare(np.ndarray)  # type: ignore
    neighbor_radius: float = cython.declare(single)  # type: ignore
    inner_radius: float = cython.declare(single)  # type: ignore

    temporal_jitter: int = cython.declare(int32_t)  # type: ignore
    rise_duration: int = cython.declare(int32_t)  # type: ignore

    decay_filtering: bool = cython.declare(bool_t)  # type: ignore
    decay_ratio: float = cython.declare(single)  # type: ignore

    localize: bool = cython.declare(bool_t)  # type: ignore

    save_shape: bool = cython.declare(bool_t)  # type: ignore
    shape_file: Optional[Path] = cython.declare(object)  # type: ignore
    cutout_start: int = cython.declare(int32_t)  # type: ignore
    cutout_end: int = cython.declare(int32_t)  # type: ignore
    cutout_length: int = cython.declare(int32_t)  # type: ignore

    verbose: bool = cython.declare(bool_t)  # type: ignore

    @cython.locals(
        recording=object,
        params=object,
        fps=single,
        l=np.ndarray,
        m=np.ndarray,
        r=np.ndarray,
        common_reference=str,
        duration_float=single,
        positions=np.ndarray,
        shape_file=object,
    )
    def __init__(self, recording: Recording, params: Params) -> None:
        self.recording = recording
        self.num_segments = recording.get_num_segments()
        self.num_frames = [
            recording.get_num_samples(seg) for seg in range(self.num_segments)
        ]
        fps = recording.get_sampling_frequency()

        self.num_channels = recording.get_num_channels()
        self.chunk_size = params["chunk_size"]

        self.rescale = params["rescale"]
        if self.rescale:
            l, m, r = np.quantile(
                self.get_random_data_chunks(), q=[0.025, 0.5, 1 - 0.025], axis=0
            )
            # quantile gives float64 on float32 data
            l: NDArray[np.single] = l.astype(np.single)
            m: NDArray[np.single] = m.astype(np.single)
            r: NDArray[np.single] = r.astype(np.single)

            self.scale: NDArray[np.single] = np.ascontiguousarray(
                params["rescale_value"] / (r - l), dtype=np.single
            )
            self.offset: NDArray[np.single] = np.ascontiguousarray(
                -m * self.scale, dtype=np.single
            )
        else:
            self.scale: NDArray[np.single] = np.ones(self.num_channels, dtype=np.single)
            self.offset: NDArray[np.single] = np.zeros(
                self.num_channels, dtype=np.single
            )

        common_reference = params["common_reference"]
        self.median_reference = common_reference == "median"
        self.average_reference = common_reference == "average"
        if self.num_channels < 20 and (self.median_reference or self.average_reference):
            warnings.warn(
                f"Number of channels too few for common {common_reference} reference"
            )

        duration_float = params["spike_duration"]
        self.spike_duration = int(duration_float * fps / 1000 + 0.5)
        duration_float = params["amp_avg_duration"]
        self.amp_avg_duration = int(duration_float * fps / 1000 + 0.5)
        self.threshold = params["threshold"]
        self.min_avg_amp = params["min_avg_amp"]
        self.max_AHP_amp = params["AHP_thr"]

        positions: NDArray[np.single] = np.array(
            [
                recording.get_channel_property(ch, "location")
                for ch in recording.get_channel_ids()
            ],
            dtype=np.single,
        )
        if positions.shape[1] > 2:
            warnings.warn(
                f"Channel locations have {positions.shape[1]} dimensions, "
                "using the last two."
            )
            positions = positions[:, -2:]
        self.positions: NDArray[np.single] = np.ascontiguousarray(
            positions, dtype=np.single
        )
        self.neighbor_radius = params["neighbor_radius"]
        self.neighbor_radius += RADIUS_EPS
        self.inner_radius = params["inner_radius"]
        self.inner_radius += RADIUS_EPS

        duration_float = params["peak_jitter"]
        self.temporal_jitter = int(duration_float * fps / 1000 + 0.5)
        duration_float = params["rise_duration"]
        self.rise_duration = int(duration_float * fps / 1000 + 0.5)

        self.decay_filtering = params["decay_filtering"]
        self.decay_ratio = params["decay_ratio"]

        self.localize = params["localize"]

        self.save_shape = params["save_shape"]
        if self.save_shape:
            shape_file = params["out_file"]
            if isinstance(shape_file, str):
                shape_file = Path(shape_file)
            shape_file.parent.mkdir(parents=True, exist_ok=True)
            if shape_file.suffix != ".bin":
                shape_file = shape_file.with_suffix(".bin")
            self.shape_file = shape_file
        else:
            self.shape_file = None
        duration_float = params["left_cutout_time"]
        self.cutout_start = int(duration_float * fps / 1000 + 0.5)
        duration_float = params["right_cutout_time"]
        self.cutout_end = int(duration_float * fps / 1000 + 0.5)
        self.cutout_length = self.cutout_start + 1 + self.cutout_end

        self.left_margin = self.temporal_jitter + max(
            self.cutout_length, self.rise_duration + 1 + self.spike_duration
        )

        self.verbose = params["verbose"]

        # sanity checks
        assert (
            self.num_channels > 0
        ), f"Expect number of channels >0, got {self.num_channels}"
        assert self.chunk_size > 0, f"Expect chunk size >0, got {self.chunk_size}"
        assert (
            self.threshold > 0
        ), f"Expect detection threshold >0, got {self.threshold}"
        assert (
            self.min_avg_amp > 0
        ), f"Expect min avg amplitude >0, got {self.min_avg_amp}"
        assert (
            self.max_AHP_amp <= 0
        ), f"Expect AHP threshold <=0, got {self.max_AHP_amp}"
        assert (
            self.neighbor_radius >= 0
        ), f"Expect neighbor radius >=0, got {self.neighbor_radius}"
        assert (
            self.inner_radius >= 0
        ), f"Expect inner neighbor radius >=0, got {self.neighbor_radius}"
        assert (
            0 <= self.decay_ratio <= 1
        ), f"Expect decay filtering ratio >=0,<=1, got {self.decay_ratio}"
        assert (
            self.temporal_jitter >= 0
        ), f"Expect temporal jitter >=0, got {self.temporal_jitter}"
        assert (
            self.spike_duration >= self.temporal_jitter
        ), f"Expect spike duration >=jitter={self.temporal_jitter}, got {self.spike_duration}"
        assert (
            self.rise_duration >= self.temporal_jitter
        ), f"Expect rising duration >=jitter={self.temporal_jitter}, got {self.rise_duration}"
        assert (
            self.cutout_start >= self.temporal_jitter
        ), f"Expect cutout start >=jitter={self.temporal_jitter}, got {self.cutout_start}"
        assert (
            self.cutout_end >= self.temporal_jitter
        ), f"Expect cutout end >=jitter={self.temporal_jitter}, got {self.cutout_end}"

    @cython.cfunc
    @cython.locals(
        chunks_per_seg=int32_t,
        chunk_size=int32_t,
        seed=object,
        chunks=list,
        seg=int32_t,
        i=int32_t,
        _random_starts=np.ndarray,
        random_starts=p_i32,
    )
    @cython.returns(np.ndarray)
    def get_random_data_chunks(
        self, chunks_per_seg: int = 20, chunk_size: int = 10000, seed: int = 0
    ) -> NDArray[np.float32]:
        # TODO: sample uniformly on samples instead of segments
        chunks: list[RealArray] = []
        for seg in range(self.num_segments):
            _random_starts = np.random.default_rng(seed).integers(  # keep a reference
                0,
                self.num_frames[seg] - chunk_size,
                size=chunks_per_seg,
                dtype=np.int32,
                endpoint=True,
            )
            random_starts = cython.cast(p_i32, _random_starts.data)
            for i in range(chunks_per_seg):
                chunks.append(
                    self.recording.get_traces(
                        segment_index=seg,
                        start_frame=random_starts[i],
                        end_frame=random_starts[i] + chunk_size,
                    )
                )
        return np.concatenate(chunks, axis=0, dtype=np.float32)

    @cython.cfunc
    @cython.locals(
        segment_index=int32_t,
        start_frame=int32_t,
        end_frame=int32_t,
        pad_left=int32_t,
        pad_right=int32_t,
        traces=np.ndarray,
        traces_float=np.ndarray,
    )
    @cython.returns(np.ndarray)
    def get_traces(
        self, segment_index: int, start_frame: int, end_frame: int
    ) -> NDArray[np.single]:
        if start_frame < 0:
            pad_left = -start_frame * self.num_channels
            start_frame = 0
        else:
            pad_left = 0
        if end_frame > self.num_frames[segment_index]:
            pad_right = (end_frame - self.num_frames[segment_index]) * self.num_channels
            end_frame = self.num_frames[segment_index]
        else:
            pad_right = 0

        traces: RealArray = self.recording.get_traces(
            segment_index=segment_index, start_frame=start_frame, end_frame=end_frame
        )

        traces_float: NDArray[np.single] = traces.astype(np.single, copy=False).reshape(
            -1
        )  # type: ignore # force cast to single

        if pad_left or pad_right:
            traces_float = np.pad(
                traces_float, (pad_left, pad_right), mode="constant", constant_values=0
            )

        return traces_float

    @cython.ccall
    @cython.returns(list)
    def detect(self) -> list[dict[str, RealArray]]:
        return [self.detect_seg(seg) for seg in range(self.num_segments)]

    @cython.cfunc
    @cython.locals(
        segment_index=int32_t,
        trace=np.ndarray,
        shape_file=object,
        num_frames=int32_t,
        chunk_start=int32_t,
        chunk_len=int32_t,
        sample_ind=np.ndarray,
        channel_ind=np.ndarray,
        amplitude=np.ndarray,
        location=np.ndarray,
        spikes=np.ndarray,
        result=dict,
    )
    @cython.returns(dict)
    def detect_seg(self, segment_index: int) -> dict[str, RealArray]:
        shape_file = (
            None
            if self.shape_file is None
            else self.shape_file.with_stem(f"{self.shape_file.stem}-{segment_index}")
        )

        det = newDet(  # type: ignore
            self.num_channels,
            self.chunk_size,
            self.left_margin,
            self.rescale,
            cython.cast(p_single, self.scale.data),
            cython.cast(p_single, self.offset.data),
            self.median_reference,
            self.average_reference,
            self.spike_duration,
            self.amp_avg_duration,
            self.threshold,
            self.min_avg_amp,
            self.max_AHP_amp,
            cython.cast(p_single, self.positions.data),
            self.neighbor_radius,
            self.inner_radius,
            self.temporal_jitter,
            self.rise_duration,
            self.decay_filtering,
            self.decay_ratio,
            self.localize,
            self.save_shape,
            str(shape_file).encode(),
            self.cutout_start,
            self.cutout_end,
        )

        num_frames = self.num_frames[segment_index]
        chunk_start = 0
        chunk_len = min(self.chunk_size, num_frames)
        while chunk_start < num_frames:
            chunk_len = min(chunk_len, num_frames - chunk_start)

            if self.verbose:
                print(
                    f"HSDetection: Analysing segment {segment_index}, "
                    f"frames from {chunk_start:8d} to {chunk_start + chunk_len:8d} "
                    f" ({100 * chunk_start / num_frames:.1f}%)"
                )

            trace = self.get_traces(
                segment_index=segment_index,
                start_frame=chunk_start - self.left_margin,
                end_frame=chunk_start + chunk_len,
            )
            det.step(cython.cast(p_single, trace.data), chunk_start, chunk_len)

            chunk_start += chunk_len

        det_len = det.finish()
        det_result = det.getResult()

        sample_ind = np.empty(det_len, dtype=np.int32)
        channel_ind = np.empty(det_len, dtype=np.int32)
        amplitude = np.empty(det_len, dtype=np.int16)
        location = np.empty((det_len, 2), dtype=np.single)
        for i in range(det_len):
            sample_ind[i] = det_result[i].frame
            channel_ind[i] = det_result[i].channel
            amplitude[i] = det_result[i].amplitude
            location[i, 0] = det_result[i].position.x
            location[i, 1] = det_result[i].position.y

        delDet(det)  # type: ignore

        if shape_file is not None and shape_file.stat().st_size > 0:
            spikes: NDArray[np.int16] = np.memmap(
                str(shape_file), dtype=np.int16, mode="r"
            ).reshape(-1, self.cutout_length)
        else:
            spikes: NDArray[np.int16] = np.empty(
                (0, self.cutout_length), dtype=np.int16
            )

        result: dict[str, RealArray] = {
            "sample_ind": sample_ind,
            "channel_ind": channel_ind,
            "amplitude": amplitude,
        }
        if self.localize:
            result["location"] = location
        if self.save_shape:
            result["spike_shape"] = spikes

        return result
