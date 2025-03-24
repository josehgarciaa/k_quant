# strategies/exact_trace.py
import numpy as np
from .trace_strategy import TraceStrategy


class RandomPhaseTrace(TraceStrategy):
    
    def __init__(self, max_iteration):
        self.max_iteration=max_iteration
    
    def __next__(self):
        """
        Returns the next trace vector. Restarts when the end is reached.

        Returns:
        -------
        np.ndarray
            Trace vector for the current iteration.
        """
        if self.iteration >= self.max_iteration:
            raise StopIteration

        x= np.random.uniform(-np.pi, np.pi, size=(self.nkpoints, self.orbdim) )
        self.iteration += 1
        print(self.iteration)
        return x