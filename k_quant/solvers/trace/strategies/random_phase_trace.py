# strategies/random_phase.py
import numpy as np
from .trace_strategy import TraceStrategy

class RandomPhaseTrace(TraceStrategy):
    def get_trace_vectors(self, D,n) -> float:
        """
        Computes the trace of the given matrix using random phase estimation.

        Args:
            matrix (np.ndarray): The input matrix.

        Returns:
            float: The estimated trace value.
        """
        random_vector = np.exp(2j * np.pi * np.random.rand(matrix.shape[0]))
        trace_estimate = np.vdot(random_vector, matrix @ random_vector)
        return np.real(trace_estimate)
