from k_quant.utils import create_mesh
import numpy as np
from .strategies.representation_strategy import SpectralRepresentationStrategy

class SpectralRepresentations:
    def __init__(self, strategy: SpectralRepresentationStrategy, spectral_operator):
        """
        Initializes the TraceComputation with a given strategy.

        Args:
            strategy (TraceStrategy): The trace computation strategy to use.
        """
        self.strategy = strategy

    def set_strategy(self, strategy: SpectralRepresentationStrategy):
        """
        Sets a new strategy for trace computation.

        Args:
            strategy (TraceStrategy): The new trace computation strategy.
        """
        self.strategy = strategy

    def get_representation_dimension(self) -> float:
        print("I am here", self.strategy)
        return self.strategy.get_representation_dimension()

    def representation_dot(self, representation_index: int, Lvector) -> float:
        print("I am here", self.strategy)
        return self.strategy.representation_dot(representation_index)


        
        

