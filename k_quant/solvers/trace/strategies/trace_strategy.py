# strategies/trace_strategy.py
from abc import ABC, abstractmethod
import numpy as np


class TraceStrategy(ABC):
    """
    Abstract base class for trace computation strategies.
    """
    nkpoints = None
    orbdim   = None
    iteration = 0
    

    def __iter__(self):
        self.iteration = 0  # Reset iteration on each new iteration start
        return self


    def reset(self):
        """
        Resets the iterator to start from the beginning.
        """
        self.iteration = 0

    def set_dimensions(self, nkpoints , orbdim):
        self.nkpoints = nkpoints
        self.orbdim = orbdim
        
    @abstractmethod
    def __next__(self):
        pass


 
