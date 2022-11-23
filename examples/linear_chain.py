from numpy import sqrt, exp, dot, conj, min,max, abs
import matplotlib.pyplot as plt
import k_quant as k

import numpy as np

print(np.show_config())


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
kdens.set_scdim(scdim=(1000,1000,1));


print("Filling the hamiltonians")
ham_k = kdens.operator( hamiltonian ); 

print("computing eigenvalues")
U_k   = kdens.eigenU(ham_k)
H_new = kdens.change_basis(U_k, ham_k);

print("diagonal elements")

eigvals= np.array([ np.diag(h) for h in H_new]).flatten();


#Chebyshev polynomials
emin = np.min(eigvals);
emax = np.max(eigvals);

alpha = 0.9;
W = (emax+emin)/2;
DE = emax - emin;
eigvals = 2*alpha*(eigvals - W )/DE; 

#DOS
M=3000;
mu = [];

print("KPM")

for m in range(M):
    T_m = np.cos( m *np.arccos(eigvals) );
    mu.append( np.sum(T_m) );

print(mu)
