# strategies/exact_trace.py
import numpy as np
from .trace_strategy import TraceStrategy

        
class ExactTrace(TraceStrategy):
    
    def __next__(self):
        """
        Returns the next trace vector. Restarts when the end is reached.

        Returns:
        -------
        np.ndarray
            Trace vector for the current iteration.
        """
        if self.iteration >= self.orbdim:
            raise StopIteration

        j = self.iteration
        x = np.zeros( (self.nkpoints, self.orbdim), dtype=np.float64)
        x[:, j::self.orbdim] = 1.0
        self.iteration += 1
        return x            
