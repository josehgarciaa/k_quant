# strategies/exact_trace.py
import numpy as np
from .representation_strategy import GreenFuncStrategy
from typing import List


class ExactRepresentation(SpectralRepresentationStrategy):
    



#    def compute(self, hamiltonian_operator: np.ndarray, energies: np.ndarray) : #returns an array of operators

#        print(compute_inverse_for_E_list(hamiltonian_operator.matrix, energies))

        
#        return 0



def compute_inverse_for_E_list(H: np.ndarray, E_list: List[np.ndarray]) -> List[np.ndarray]:
    """
    Computes the inverse of (H - E) for each E in the input list,
    where H and each E are block-diagonal matrices with shape (D, n, n).
    
    For each E in E_list, this function computes:
        (H - E)^(-1)
    by subtracting E from H blockwise and then inverting each resulting block-diagonal matrix.
    
    Parameters:
        H (np.ndarray): A block-diagonal matrix of shape (D, n, n).
        E_list (List[np.ndarray]): A list of block-diagonal matrices, each of shape (D, n, n).
        
    Returns:
        List[np.ndarray]: A list of block-diagonal matrices, each of shape (D, n, n),
                          where each element is the inverse of (H - E) for the corresponding E.
    """
    inverses = []
    for idx, E in enumerate(E_list):
        diff = H - E  # diff is shape (D, n, n)
        try:
            # Try to use batched inversion if available
            inv_diff = np.linalg.inv(diff)
        except np.linalg.LinAlgError:
            # Fall back to inverting each block individually
            D, n, _ = diff.shape
            inv_diff = np.array([np.linalg.inv(diff[d]) for d in range(D)])
        inverses.append(inv_diff)
    return inverses

class ExactRepresentation(GreenFuncStrategy):
    
    def compute(self, hamiltonian_operator: np.ndarray, energies: np.ndarray) : #returns an array of operators

        print(compute_inverse_for_E_list(hamiltonian_operator.matrix, energies))

        
        return 0
