import numpy as np
from k_quant.utils.ham_function_optimizer import optimize_k_hamiltonian, compute_hamiltonian
from k_quant.global_parameters import get_kmesh

class Operator:

    def __init__(self, op_function, name="generic_operator"):
        self.name= name
        self.op_function =op_function
        self.create_operator()
        self.shape = self.matrix.shape

    def GetMatrix(self):
        return self.matrix

    def create_operator(self):
        get_kmesh()
        self.matrix = compute_hamiltonian(get_kmesh(), self.op_function)
        return self

    def dot(self, other):
        self.op_function(other)
        
        
        
                  

