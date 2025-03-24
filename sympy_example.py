import matplotlib.pyplot as plt
#import k_quant as k

from k_quant.global_parameters import set_kmesh, get_kmesh, set_energy_grid, get_energy_grid
import numpy as np
from k_quant.operators.operator import Operator
from k_quant.operators.spectral_operator import SpectralOperator

from k_quant.solvers.trace.trace_calculation import Trace 
from k_quant.solvers.trace.strategies import  ExactTrace

from k_quant.solvers.spectral_representations.spectral_representation import SpectralRepresentations
from k_quant.solvers.spectral_representations.strategies.chebyshev import ChebyshevRepresentation 

from k_quant.operators.spectral_operators_strategies.spectra_operator_strategy import SpectralOperatorStrategy
from k_quant.operators.spectral_operators_strategies.advanced_green_function import AdvancedGreenFuntion
from k_quant.operators.spectral_operators_strategies.retarded_green_function import RetardedGreenFuntion
from k_quant.operators.spectral_operators_strategies.derivative_advanced_green_function import DerivateAdvancedGreenFuntion
from k_quant.operators.spectral_operators_strategies.derivative_retarded_green_function import DerivateRetardedGreenFuntion
from k_quant.operators.spectral_operators_strategies.ImGreenFunction import ImAdvancedGreenFunction

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("tkagg")  # Use TkAgg backend
plt.ion()  # Turn on interactive mode


"""
Let us first define a Hamiltonian in momentum space
and for this purpose we choose as prototype a model for p_z electrons
 in graphene within the nearest neighbor approximation with lattice vectors
 defined in the lat_vec variable
"""
lat_vec = [ [ 1, 0,0 ], [ 0,1,0 ], [ 0,0,1 ] ]

def Ham_k(k): #MUST BE IN RECIPROCAL
    a_0, a_1, a2  = lat_vec
    hop = 1.0
    f_k = 2*hop* np.cos( 2*np.pi*k[0]*a_0[0] ) #
    return np.array([ [ f_k] ])

def Vel_k(k): #MUST BE IN RECIPROCAL
    a_0, a_1, a2  = lat_vec
    hop = 1.0
    f_k = 4*hop*a_0[0]* np.cos( 2*np.pi*k[0]*a_0[0] ) #
    return np.array([ [ f_k] ])





import sympy as sp
import numpy as np
from numba import njit

# 1. Define the symbolic Hamiltonian in 2D k-space.
kx, ky, t = sp.symbols('kx ky t', real=True)
I = sp.I  # imaginary unit

# Define the off-diagonal element: t*(exp(i*kx) + exp(i*ky))
f = t * (sp.exp(I * kx) + sp.exp(I * ky))

# Construct the 2x2 Hamiltonian
H_sym = sp.Matrix([[0, f],
                   [sp.conjugate(f), 0]])

# Compute the symbolic derivatives.
dH_dkx_sym = H_sym.diff(kx)
dH_dky_sym = H_sym.diff(ky)

# Lambdify the Hamiltonian and its derivatives for numerical evaluation.
H_func     = sp.lambdify((kx, ky, t), H_sym, modules='numpy')
dH_dkx_func = sp.lambdify((kx, ky, t), dH_dkx_sym, modules='numpy')
dH_dky_func = sp.lambdify((kx, ky, t), dH_dky_sym, modules='numpy')

# Test the lambdified functions for a sample kx, ky, and t.
kx_val = 0.5
ky_val = 0.3
t_val  = 1.0




print("Hamiltonian using lambdify at kx=0.5, ky=0.3, t=1.0:")
print(H_func(kx_val, ky_val, t_val))
print("\nDerivative with respect to kx (lambdified):")
print(dH_dkx_func(kx_val, ky_val, t_val))
print("\nDerivative with respect to ky (lambdified):")
print(dH_dky_func(kx_val, ky_val, t_val))

import inspect
source_code = inspect.getsource(H_func)


energies = np.linspace(-3,3,100)
set_energy_grid(energies)
set_kmesh( (10000, 1 , 1))


