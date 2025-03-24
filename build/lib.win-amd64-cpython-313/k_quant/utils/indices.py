import numpy as np

def create_iteration_indices(sp_ops_repr_dim):
    """
    Create iteration indices for nested loops based on dimensions provided.

    Parameters:
    -----------
    sp_ops_repr_dim: list or array-like
        A list containing dimensions for each loop.

    Returns:
    --------
    indices: ndarray
        An array of shape (total_combinations, number_of_loops) representing indices.
    grids: list of ndarray
        List containing meshgrid arrays corresponding to each dimension.
    """
    ranges = [np.arange(d) for d in sp_ops_repr_dim]
    grids = np.meshgrid(*ranges, indexing='ij')
    indices = np.array([grid.flatten() for grid in grids]).T

    return indices