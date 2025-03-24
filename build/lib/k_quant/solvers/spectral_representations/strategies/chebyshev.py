# strategies/exact_trace.py
import numpy as np
from .representation_strategy import SpectralRepresentationStrategy
from typing import List


class ChebyshevRepresentation(SpectralRepresentationStrategy):

    def get_representation_dimension(self) -> int:
        print("chebyshev representation computed using eta")
        return 3

    def representation_dot(self, representation_index: int, Lvector) -> float:
        print("compute the Chebyshev Vector for", representation_index)
        return Lvector
