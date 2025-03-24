import numpy as np
from k_quant.operators.operator import Operator
from k_quant.operators.spectral_operators_strategies.spectra_operator_strategy import SpectralOperatorStrategy
from k_quant.operators.spectral_operators_strategies.advanced_green_function import AdvancedGreenFuntion
from k_quant.operators.spectral_operators_strategies.retarded_green_function import RetardedGreenFuntion
from k_quant.operators.spectral_operators_strategies.derivative_advanced_green_function import DerivateAdvancedGreenFuntion
from k_quant.operators.spectral_operators_strategies.derivative_retarded_green_function import DerivateRetardedGreenFuntion
from k_quant.operators.spectral_operators_strategies.ImGreenFunction import ImGreenFunction



class SpectralOperator:

    def __init__(self,  strategy: SpectralOperatorStrategy,
                        hamiltonian_op, 
                        broadening, 
                        name="spectral_operator"):
        self.name= name
        self.strategy = strategy        
        self.hamiltonian_op = hamiltonian_op
        self.broadening = broadening
        self.shape = hamiltonian_op.matrix.shape


    def set_strategy(self, strategy: SpectralOperatorStrategy):
        """
        Sets a new strategy for trace computation.

        Args:
            strategy (TraceStrategy): The new trace computation strategy.
        """
        self.strategy = strategy
        
    def GetMatrix(self, energy):
        return self.strategy.GetMatrix(self.hamiltonian_op, self.broadening,energy)

        
