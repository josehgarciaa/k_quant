import numpy as np
from k_quant.operators.spectral_operators_strategies.spectra_operator_strategy import SpectralOperatorStrategy

class RetardedGreenFuntion(SpectralOperatorStrategy):

    def GetMatrix(self, hamiltonian_op, broadening, energy):
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
        hamiltonian_matrix = np.array(hamiltonian_op.GetMatrix())  # Shape (D, n, n)
        n = hamiltonian_matrix.shape[1]
                
        I = np.eye(n)[np.newaxis, :, :]
        shifted_H = I*(energy + 1j * broadening) -hamiltonian_matrix    # Shape: (D, n, n)
        green_matrix =-np.linalg.inv(shifted_H )

        return green_matrix


        
    def SpectralFunction(self, x, broadening, energy):
        return 1/(energy -x + 1j*broadening)