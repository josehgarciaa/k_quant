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




energies = np.linspace(-3,3,100)
set_energy_grid(energies)
set_kmesh( (10000, 1 , 1))


hamiltonian_op = Operator(op_function=Ham_k, name="Hamiltonian" )
velocity_op = Operator(op_function = Vel_k, name="VX" )


eta = 0.001
green_func_op = SpectralOperator( strategy = ImAdvancedGreenFunction(),
                                  hamiltonian_op =hamiltonian_op,  
                                  broadening = eta,
                                  name="AdvancedGreen's Function" )


my_trace  = Trace(ExactTrace())
#cheb_green_func_op= SpectralRepresentations(ChebyshevRepresentation(), green_func_op)
#dos_theo = np.array([ 1/np.sqrt(4- E*E ) for E in get_energy_grid() ])*10000
#dos = my_trace.compute( green_func_op)

#plt.plot(get_energy_grid(), dos)
#plt.plot(get_energy_grid(), dos_theo)

#plt.show()
#input("Press Enter to exit...")  # Keeps the window open


my_trace  = Trace(ExactTrace())
cond = my_trace.compute(  velocity_op, green_func_op, velocity_op, green_func_op)

plt.plot(get_energy_grid(), cond)

plt.show()
input("Press Enter to exit...")  # Keeps the window open
