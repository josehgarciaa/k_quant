import numpy as np


def batch_dot_product(matrices, vectors):
    """
    Efficient dot product between a collection of matrices and vectors.

    Parameters:
    -----------
    matrices : np.ndarray
        Array of shape (D, n, n) containing D matrices of size n x n.
    vectors : np.ndarray
        Array of shape (D, n) containing D vectors of size n.

    Returns:
    --------
    np.ndarray
        An array of shape (D, n) containing the dot products.
    """
    
    
    return np.array([ H.dot(x) for H,x in zip(matrices, vectors)])



def vectorized_apply_A(blocks: np.ndarray, x: np.ndarray) -> np.ndarray:
    """
    Applies a block-diagonal operator (given as an array of blocks)
    to a vector or matrix x in a vectorized fashion.
    
    Parameters:
        blocks (np.ndarray): Array of shape (D, n, n) where each blocks[d] is an n×n matrix.
        x (np.ndarray): A vector of shape (D*n,) or a matrix of shape (D*n, k).
        
    Returns:
        np.ndarray: The result of applying the operator to x, with the same shape as x.
    """
    D, n, _ = blocks.shape
    if x.ndim == 1:
        # Reshape x to (D, n, 1)
        x_block = x.reshape(D, n, 1)
        y_block = np.matmul(blocks, x_block)  # (D, n, 1)
        return y_block.reshape(D * n)
    elif x.ndim == 2:
        # Assume shape is (D*n, k)
        k = x.shape[1]
        x_block = x.reshape(D, n, k)
        y_block = np.matmul(blocks, x_block)  # (D, n, k)
        return y_block.reshape(D * n, k)
    else:
        raise ValueError("x must be a 1D or 2D array.")
