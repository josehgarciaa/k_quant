# strategies/exact_trace.py
import numpy as np
from .trace_strategy import TraceStrategy


class TraceVectorGenerator:
    def __init__(self, D, n):
        """
        Initializes the trace vector generator.

        Parameters:
        -----------
        D: int
            Dimension size.
        n: int
            Number of trace vectors to generate.
        """
        self.D = D
        self.n = n
        self.iteration = 0

    def __iter__(self):
        self.iteration = 0  # Reset iteration on each new iteration start
        return self

    def __next__(self):
        """
        Returns the next trace vector. Restarts when the end is reached.

        Returns:
        -------
        np.ndarray
            Trace vector for the current iteration.
        """
        if self.iteration >= self.n:
            raise StopIteration

        j = self.iteration
        x = np.zeros( (self.D, self.n), dtype=np.float64)
        x[:, j::self.n] = 1.0
        self.iteration += 1
        return x

    def reset(self):
        """
        Resets the iterator to start from the beginning.
        """
        self.iteration = 0



class ExactTrace(TraceStrategy):
    
    def get_trace_vectors(self, D, n):
        """
        Generator to yield trace vectors for each iteration.

        Parameters:
        -----------
        D: int
            Dimension size.
        n: int
            Number of trace vectors to generate.

        Yields:
        -------
        np.ndarray
            Trace vector for each iteration.
            
        """
        
        generator = TraceVectorGenerator(D,n)
        return generator
            
