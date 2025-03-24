from k_quant.utils import create_mesh
import numpy as np
from .strategies.trace_strategy import TraceStrategy

from k_quant.operators.operator import Operator
from k_quant.operators.spectral_operator import SpectralOperator
from k_quant.solvers.spectral_representations.spectral_representation import SpectralRepresentations
from .compute_methods import exact_compute, representation_compute
from k_quant.utils.validators import validate_operator_dimensions



class Trace:
    def __init__(self, strategy: TraceStrategy):
        """
        Initializes the TraceComputation with a given strategy.

        Args:
            strategy (TraceStrategy): The trace computation strategy to use.
        """
        self.strategy = strategy

    def set_strategy(self, strategy: TraceStrategy):
        """
        Sets a new strategy for trace computation.

        Args:
            strategy (TraceStrategy): The new trace computation strategy.
        """
        self.strategy = strategy


    def get_trace_vectors(self, nkpoints, orbdim) -> float:
        """
        Computes the trace using the current strategy.

        Args:
            matrix (np.ndarray): The input matrix.

        Returns:
            float: The computed trace.
        """
        print("I got the trace vector for my strategy")
        self.strategy.set_dimensions(nkpoints, orbdim)
        return self.strategy


    def compute(self, *operators) -> float:
        """
        Computes the trace using the current strategy.

        Args:
            matrix (np.ndarray): The input matrix.

        Returns:
            float: The computed trace.
        """


        

        #Here we identify if the spectral are in a representation or exact
        is_exact = False
        is_sp_representation = False
        illegal_type =False
        
        operator_list = []
        for op in operators:
            if isinstance(op, Operator) or isinstance(op, SpectralOperator):
                operator_list.append(op)
                is_exact =True
            elif isinstance(op, SpectralRepresentations):
                operator_list.append(op)
                is_sp_representation =True
                is_exact =False
            else:
                illegal_type =True


        D, n, _ = validate_operator_dimensions( operators )
        trace_vectors = self.get_trace_vectors(D,n) 


        if is_exact:
            return exact_compute(operator_list,trace_vectors)

        if is_sp_representation:
            return representation_compute(operator_list, trace_vectors)

        if illegal_type:
            raise("Illegal Object passed to compute function in Trace")


        #trace_vector = self.get_trace_vectors(self, matrix: np.ndarray)
        
        #return self.strategy.compute_trace(matrix)
    