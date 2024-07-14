from typing import (Any, Iterable, Literal, Optional, Protocol, Union,
                    runtime_checkable)

import numpy as np
from numpy.typing import NDArray

__all__ = ['Recording', 'RealArray']


RealArray = Union[NDArray[np.integer], NDArray[np.floating]]


@runtime_checkable
class Recording(Protocol):
    """Interface for recording data, of which spikeinterface.core.BaseRecording
    is a virtual subclass without further registration.
    """

    def get_num_channels(self) -> int: ...

    def get_channel_ids(self) -> Iterable[Any]: ...

    def get_channel_property(self, channel_id: Any, key: Any) -> Any: ...

    def get_sampling_frequency(self) -> float: ...

    def get_num_segments(self) -> int: ...

    def get_num_samples(self, segment_index: Optional[int] = None) -> int: ...

    def get_traces(self,
                   segment_index: Optional[int] = None,
                   start_frame: Optional[int] = None,
                   end_frame: Optional[int] = None,
                   channel_ids: Optional[Iterable[Any]] = None,
                   order: Optional[Literal['C', 'F']] = None,
                   return_scaled: bool = False
                   ) -> RealArray: ...
