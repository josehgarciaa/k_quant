from numpy import sqrt, exp, dot, conj, min,max, abs
import matplotlib.pyplot as plt
import k_quant as k
from k_quant import operators as ops

import numpy as np
"""
Let us first define a Hamiltonian in momentum space
and for this purpose we choose as prototype a model for p_z electrons
 in graphene within the nearest neighbor approximation with lattice vectors
 defined in the lat_vec variable
"""
lat_vec = [ [ 1/2, sqrt(3)/2,0 ], [ 1/2,-sqrt(3)/2,0 ], [ 0,0,1 ] ]
def hamiltonian(k):
    a_0, a_1, a2  = lat_vec;
    hop = 2.8;
    f_k = hop*( 1 + exp( -1j*dot(k,a_0)) + exp( -1j*dot(k,a_1)) );
    return [ [ 0        , f_k],
             [ conj(f_k),  0 ]
            ];

#The band class requires lattice vectors and the hamiltonian function
kdens = k.Density(lat_vec, hamiltonian );
kdens.set_scdim(scdim=(10,10,1));

ops.energy
ham_k = kdens.operator( hamiltonian ); 

print( np.eigvals.shape );
