# strategies/trace_strategy.py
from abc import ABC, abstractmethod
import numpy as np

class TraceStrategy(ABC):
    """
    Abstract base class for trace computation strategies.
    """

    @abstractmethod
    def get_trace_vectors(self, D,n) -> float:
        """
        Computes the trace of the given matrix.

        Args:
            matrix (np.ndarray): The input matrix.

        Returns:
            float: The computed trace value.
        """
        pass
