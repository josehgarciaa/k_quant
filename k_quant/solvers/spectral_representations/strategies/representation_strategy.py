# strategies/trace_strategy.py
from abc import ABC, abstractmethod
import numpy as np

class SpectralRepresentationStrategy(ABC):
    """
    Abstract base class for trace computation strategies.
    """

    @abstractmethod
    def get_representation_dimension(self) -> float:
        pass

    @abstractmethod
    def representation_dot(self, representation_index: int, Lvector ) -> float:
        pass

