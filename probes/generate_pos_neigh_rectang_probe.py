from itertools import product

import numpy as np
from scipy.spatial.distance import cdist

name_of_probe = "neuroseeker_128"
n_channels = 128
n_rows = 32
n_cols = 4

assert n_rows * n_cols == n_channels, "Product of the number of rows and "\
    "columns doesn't equal number of channels."

ch_positions = np.asarray(list(product(range(n_rows), range(n_cols))))

# NB: Notice the column, row order in write
with open("positions_{}".format(name_of_probe), 'w') as f:
    for pos in ch_positions:
        f.write('{},{},\n'.format(pos[1], pos[0]))

# NB: it is also possible to use metric='cityblock' (Manhattan distance)
distances = cdist(ch_positions, ch_positions, metric='euclidean')
radius = 2.9

indices = np.arange(n_channels)
with open("neighbormatrix_{}".format(name_of_probe), "w") as f:
    for dist_from_ch in distances:
        neighbors = indices[dist_from_ch <= radius]
        f.write('{},\n'.format(str(list(neighbors))[1:-1]))
