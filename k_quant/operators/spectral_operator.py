import numpy as np
from k_quant.operators.operator import Operator

class SpectralOperator:

    def __init__(self,  hamiltonian_op, 
                        broadening, 
                        spectral_function, 
                        name="spectral_operator"):
        self.name= name
        self.hamiltonian_op = hamiltonian_op
        self.broadening = broadening
        self.spectral_function = spectral_function 
        self.shape = hamiltonian_op.matrix.shape

    def dot(self, other):
        print("implement spectral doct")

    def GetMatrix(self, energy):
        """
        Efficiently calculates the Green's function matrix using broadcasting.

        Parameters:
        -----------
        energy : float
            The energy value used to compute the matrix.

        Returns:
        --------
        np.ndarray
            The Green's function matrix with shape (D, n, n).
        """
        hamiltonian_matrix = np.array(self.hamiltonian_op.GetMatrix())  # Shape (D, n, n)
        n = hamiltonian_matrix.shape[1]
                
        # Create the diagonal adjustment using broadcasting
        diag_adjustment = np.eye(n) * (energy + 1j * self.broadening)

        # Use broadcasting to subtract the diagonal adjustment from each block
        green_matrix = hamiltonian_matrix - diag_adjustment[np.newaxis, :, :]
        green_matrix = np.array([ np.linalg.inv(G) for G in green_matrix ])     #OPTIMIZATION PROBLEM
        

        return green_matrix

        
