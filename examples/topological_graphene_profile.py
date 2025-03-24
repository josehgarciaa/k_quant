import cProfile
import pstats
import io
import numpy as np

# Import the necessary k_quant modules
import k_quant.global_parameters as kparam
import k_quant.operators as kop  
import k_quant.solvers.trace as ktr  
import k_quant.solvers.trace.strategies as ktr_strategy  
import k_quant.operators.spectral_operators_strategies as ksp_type
from k_quant.utils.k_operator_optimizer import optimize_k_operator
from k_quant.utils.calculus import cumulative_integral

# Define the lattice vectors and parameters
lat_vec = np.array([[3/2, np.sqrt(3)/2, 0],
                    [3/2,-np.sqrt(3)/2, 0],
                    [0,0,1]])

t = 2.8   # Nearest-neighbor hopping energy (eV)
lambda_soc = 1.0/3/np.sqrt(3)  # Next-nearest-neighbor hopping energy (eV)

# Nearest-neighbor vectors
delta1 = np.array([1/2, np.sqrt(3)/2, 0])
delta2 = np.array([1/2,-np.sqrt(3)/2, 0])
delta3 = np.array([-1, 0, 0])
delta = [delta1, delta2, delta3]

# Next-nearest-neighbor vectors
b1 = lat_vec[0]
b2 = lat_vec[1] - lat_vec[0]
b3 = lat_vec[1]

# Hamiltonian function
def Ham_k(k):
    H = np.zeros((2, 2), dtype=complex)
    # Nearest-neighbor contribution
    gamma = sum(np.exp(1j * (k.dot(d))) for d in delta)
    H[0, 1] = -t * gamma
    H[1, 0] = np.conjugate(H[0, 1])
    # Next-nearest-neighbor contribution
    gamma_prime = 2*(np.sin(k.dot(b1)) + np.sin(k.dot(b2)) - np.sin(k.dot(b3)))
    H[0, 0] = -lambda_soc * gamma_prime
    H[1, 1] = -H[0, 0]
    return H

# Velocity in x-direction
def Vel_x(k):
    H = np.zeros((2, 2), dtype=complex)
    gamma = 1j * sum(d[0]*np.exp(1j * (k.dot(d))) for d in delta)
    H[0, 1] = -t * gamma
    H[1, 0] = np.conjugate(H[0, 1])
    gamma_prime = 2*(b1[0]*np.cos(k.dot(b1)) + b2[0]*np.cos(k.dot(b2)) - b3[0]*np.cos(k.dot(b3)))
    H[0, 0] = -lambda_soc * gamma_prime
    H[1, 1] = -H[0, 0]
    return H

# Velocity in y-direction
def Vel_y(k):
    H = np.zeros((2, 2), dtype=complex)
    gamma = 1j * sum(d[1]*np.exp(1j * (k.dot(d))) for d in delta)
    H[0, 1] = -t * gamma
    H[1, 0] = np.conjugate(H[0, 1])
    gamma_prime = 2*(b1[1]*np.cos(k.dot(b1)) + b2[1]*np.cos(k.dot(b2)) - b3[1]*np.cos(k.dot(b3)))
    H[0, 0] = -lambda_soc * gamma_prime
    H[1, 1] = -H[0, 0]
    return H

def main():
    # Set up the simulation parameters
    n0, n1 = 100, 100 
    energies = np.linspace(-13, 13, 1000)
    kparam.set_lattice_vector(lat_vec)
    kparam.set_energy_grid(energies)
    kparam.set_kmesh((n0, n1, 1))
    
    # Create operators
    hamiltonian_op = kop.Operator(op_function=Ham_k, name="Hamiltonian")
    velX_op = kop.Operator(op_function=Vel_x, name="VelocityX")
    velY_op = kop.Operator(op_function=Vel_y, name="VelocityY")
    
    broadening = 0.1
    ImGF_op = kop.SpectralOperator(strategy=ksp_type.ImGreenFunction(),
                                   hamiltonian_op=hamiltonian_op,  
                                   broadening=broadening,
                                   name="ImG's Function")
    
    adv_DGF_op = kop.SpectralOperator(strategy=ksp_type.DerivateAdvancedGreenFuntion(),
                                      hamiltonian_op=hamiltonian_op,  
                                      broadening=broadening,
                                      name="Derivative AdvancedGreen's Function")
    
    my_trace = ktr.Trace(ktr_strategy.ExactTrace())
    
    # Compute trace for different operator combinations
    condxx_KERNEL = my_trace.compute(velX_op, adv_DGF_op, velX_op, ImGF_op)
    condxy_KERNEL = my_trace.compute(velX_op, adv_DGF_op, velY_op, ImGF_op)
    
    # For demonstration, we'll print summaries of the computed results.
    print("condxx_KERNEL summary:")
    print(condxx_KERNEL)
    
    print("\ncondxy_KERNEL summary:")
    print(condxy_KERNEL)

if __name__ == "__main__":
    profiler = cProfile.Profile()
    profiler.enable()
    
    main()  # Execute the main function
    
    profiler.disable()
    
    # Create a stream to capture profiling statistics
    stream = io.StringIO()
    stats = pstats.Stats(profiler, stream=stream).sort_stats("cumtime")
    stats.print_stats()
    
    # Output profiling results to the console
    print("\nProfiling results:")
    print(stream.getvalue())
