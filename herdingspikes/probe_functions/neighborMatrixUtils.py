import numpy as np


def createNeighborMatrix(neighbor_matrix_name, positions_file_path, neighbor_radius):
    position_file = open(positions_file_path, "r")
    positions = []
    for position in position_file.readlines():
        positions.append(np.array(position[:-2].split(",")).astype(float))
    channel_positions = np.asarray(positions)
    NUM_CHANNELS = len(channel_positions)

    neighbor_matrix = []
    for channel in range(NUM_CHANNELS):
        # Calculate distance from current channel to all other channels
        curr_channel_distances = np.sqrt(
            np.sum((channel_positions - channel_positions[channel]) ** 2, axis=1)
        )

        # Find all neighbors in given radius and add them to neighbors
        neighbors = np.where(curr_channel_distances < neighbor_radius)[0]
        neighbor_matrix.append(neighbors)
    position_file.close()
    writeoutNeighborMatrix(neighbor_matrix, neighbor_matrix_name)


def writeoutNeighborMatrix(neighbor_matrix, neighbor_matrix_name):
    neighbor_matrix_file = open(neighbor_matrix_name, "w")
    for channel_number, neighbors in enumerate(neighbor_matrix):
        for neighbor in neighbors:
            neighbor_matrix_file.write(str(neighbor))
            neighbor_matrix_file.write(str(","))
        neighbor_matrix_file.write("\n")
    neighbor_matrix_file.close()
