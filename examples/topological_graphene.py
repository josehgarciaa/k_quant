import matplotlib.pyplot as plt
import numpy as np
import numpy as np

import k_quant.global_parameters as kparam
import k_quant.operators as kop  
import k_quant.solvers.trace as ktr  
import k_quant.solvers.trace.strategies as ktr_strategy  

import k_quant.operators.spectral_operators_strategies as ksp_type
from k_quant.utils.k_operator_optimizer import  optimize_k_operator


import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("tkagg")  # Use TkAgg backend
plt.ion()  # Turn on interactive mode


import k_quant.global_parameters as kparam

from k_quant.utils.calculus import cumulative_integral

lat_vec = np.array([ [3/2, np.sqrt(3)/2, 0], [3/2,-np.sqrt(3)/2, 0], [0,0,1]])

t = 2.8   # Nearest-neighbor hopping energy (eV)
lambda_soc = 1.0/3/np.sqrt(3)  # Next-nearest-neighbor hopping energy (eV)
                                            
    # Nearest-neighbor vectors
delta1 = np.array([1/2, np.sqrt(3)/2, 0])
delta2 = np.array([1/2,-np.sqrt(3)/2, 0])
delta3 = np.array([-1, 0, 0])

delta = [delta1, delta2, delta3]

    # Next-nearest-neighbor vectors
b1 = lat_vec[0]
b2 = lat_vec[1]-lat_vec[0]
b3 = lat_vec[1]


# Hamiltonian function
def Ham_k(k):

    H = np.zeros((2, 2), dtype=complex)

    # Nearest-neighbor contribution
    gamma = sum(np.exp( 1j * (k.dot(d) )) for d in delta)
    H[0, 1] = -t * gamma
    H[1, 0] = np.conjugate(H[0, 1])

    # Next-nearest-neighbor contribution
    gamma_prime = 2*( np.sin(k.dot(b1))+ np.sin(k.dot(b2)) - np.sin(k.dot(b3)) )
    H[0, 0] = -lambda_soc * gamma_prime
    H[1, 1] = -H[0, 0]

    return H


# Hamiltonian function
def Vel_x(k):

    H = np.zeros((2, 2), dtype=complex)

    # Nearest-neighbor contribution
    gamma = 1j*sum( d[0]*np.exp( 1j * (k.dot(d) )) for d in delta)
    H[0, 1] = -t * gamma
    H[1, 0] = np.conjugate(H[0, 1])

    # Next-nearest-neighbor contribution
    gamma_prime = 2*( b1[0]*np.cos(k.dot(b1))+ b2[0]*np.cos(k.dot(b2)) - b3[0]*np.cos(k.dot(b3)) )
    H[0, 0] = -lambda_soc * gamma_prime
    H[1, 1] = -H[0, 0]

    return H


# Hamiltonian function
def Vel_y(k):

    H = np.zeros((2, 2), dtype=complex)

    # Nearest-neighbor contribution
    gamma = 1j*sum( d[1]*np.exp( 1j * (k.dot(d) )) for d in delta)
    H[0, 1] = -t * gamma
    H[1, 0] = np.conjugate(H[0, 1])

    # Next-nearest-neighbor contribution
    gamma_prime = 2*( b1[1]*np.cos(k.dot(b1))+ b2[1]*np.cos(k.dot(b2)) - b3[1]*np.cos(k.dot(b3)) )
    H[0, 0] = -lambda_soc * gamma_prime
    H[1, 1] = -H[0, 0]

    return H


n0,n1 =300, 300 
energies = np.linspace(-13,13,1000)
kparam.set_lattice_vector(lat_vec)
kparam.set_energy_grid(energies)
kparam.set_kmesh( (n0, n1 , 1))


hamiltonian_op = kop.Operator(op_function=Ham_k, name="Hamiltonian" )
velX_op = kop.Operator(op_function=Vel_x, name="VelocityX" )
velY_op = kop.Operator(op_function=Vel_y, name="VelocityY" )


broadening = 0.1
ImGF_op = kop.SpectralOperator( strategy = ksp_type.ImGreenFunction(),
                                hamiltonian_op =hamiltonian_op,  
                                broadening = broadening,
                                name="ImG's Function" )

adv_DGF_op = kop.SpectralOperator(  strategy = ksp_type.DerivateAdvancedGreenFuntion(),
                                    hamiltonian_op =hamiltonian_op,  
                                    broadening = broadening,
                                    name="Derivative AdvancedGreen's Function" )





my_trace  = ktr.Trace(ktr_strategy.ExactTrace())
energies = kparam.get_energy_grid()

#condxx_KERNEL = my_trace.compute(  velX_op, adv_DGF_op, velX_op, ImGF_op)
#plt.plot(energies, np.pi*cumulative_integral(energies, np.imag(condxx_KERNEL /n0/n1)), label='case_1')

#condxx = my_trace.compute(  velX_op, ImGF_op, velX_op, ImGF_op)
#plt.plot(energies,  np.real(condxx) /n0/n1, label='case_1')


condxy_KERNEL = my_trace.compute(  velX_op, adv_DGF_op, velY_op, ImGF_op)
#plt.plot(energies, np.pi*cumulative_integral(energies, np.imag(condxy_KERNEL /n0/n1)), label='case_1')


plt.show()
input("Press Enter to exit...")  # Keeps the window open

