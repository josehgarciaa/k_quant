
from k_quant.global_parameters import get_energy_grid

from k_quant.operators.spectral_operator import SpectralOperator
from k_quant.operators.operator import Operator
from k_quant.utils.indices import create_iteration_indices
from k_quant.solvers.spectral_representations.spectral_representation import SpectralRepresentations 

from k_quant.linalg.vectorized_mat_vec import batch_dot_product, vectorized_apply_A
from k_quant.utils.validators import validate_operator_dimensions
import copy
import numpy as np 

def representation_compute(operators,trace_vectors):
    print("Representation compute")


def _representation_compute(operators,trace_vectors):

    operators_list = []
    sp_ops_repr_dim = []
    # Loop over arguments to gather expansions
    for op in operators:
        if isinstance(op, SpectralOperator):
            sp_rep_op = SpectralRepresentations(op)
            operators_list.append(sp_rep_op)
            sp_ops_repr_dim.append( sp_rep_op.get_representation_dimension())
        elif isinstance(op, Operator) or isinstance(op, Operator):
            operators_list.append(op)
        else:
            raise TypeError("Arguments must be Operators or SpectralOperators")

    entire_repr_indices = create_iteration_indices( sp_ops_repr_dim )
    representation_matrix = np.zeros(entire_repr_indices.shape, dtype=complex)

    for trace_vector in trace_vectors:
        temporal_vector = trace_vector 
        for repr_indices in entire_repr_indices:

            for op_idx, op in enumerate(operators_list):
                if isinstance(op, SpectralRepresentations):
                    temporal_vector = op.representation_dot(repr_indices[op_idx],temporal_vector)
                elif isinstance(arg, Operator):
                    temporal_vector = op.dot(temporal_vector)

            representation_matrix[*repr_indices]= trace_vector.dot(temporal_vector)        
    
    

def exact_compute(operators,trace_vectors):
    energies = get_energy_grid()
 
    Tr_Es = np.zeros(energies.shape, dtype=complex)  
    for e_idx, energy in enumerate(energies):
        #print(e_idx, energy)
        Tr_Es[e_idx] = 0.0
        trace_vectors.reset()
        for x in trace_vectors:
            y = copy.copy(x)
            for op in operators:
                if isinstance(op, Operator):
                    this_matrix = op.GetMatrix()
                elif isinstance(op, SpectralOperator):
                    this_matrix = op.GetMatrix(energy)
                else:
                    raise ValueError(f"Unsupported operator type: {type(op)}")
                
                # Use vectorized_apply_A.
                # Since vectorized_apply_A expects x to be either a flat array (if 1D)
                # or a (D*n, k) 2D array, we first reshape y to a flat array.
                D, n = y.shape
                y_flat = y.reshape(-1)  # Shape (D*n,)
                # Apply the block-diagonal operator in a vectorized way.
                y_flat = vectorized_apply_A(this_matrix, y_flat)
                # Reshape the result back to (D, n)
                y = y_flat.reshape(D, n)

#                    y = batch_dot_product(op.GetMatrix(energy), y)                    
            Tr_Es[e_idx] += np.vdot(x.flatten(), y.flatten())
    return Tr_Es
 
 

