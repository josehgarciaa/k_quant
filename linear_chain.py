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
    return np.array([ [ f_k, 0], [ 0, f_k ] ])
    

    
def green_function(x, E, eta):
    return 1.0/(x - E + 1j*eta)


energies = np.linspace(-3,3,100)
set_energy_grid(energies)
set_kmesh( (10000, 1 , 1))


hamiltonian_op = Operator(op_function=Ham_k, name="Hamiltonian" )

print(hamiltonian_op.GetMatrix())

eta = 0.001
green_func_op = SpectralOperator( hamiltonian_op =hamiltonian_op,  
                                  broadening = eta,
                                  spectral_function = green_function, name="AdvancedGreen's Function" )


my_trace  = Trace(ExactTrace())
cheb_green_func_op= SpectralRepresentations(ChebyshevRepresentation(), green_func_op)

print( [x[0] for x in get_kmesh()])

my_cos = [np.real(H[0,0]) for H in hamiltonian_op.GetMatrix()]
my_invcos_exact = [eta/(H[0,0]*H[0,0] + eta*eta) for H in hamiltonian_op.GetMatrix()]

dos_exact = np.array([ np.sum([np.imag(1/(H[0,0] -(E + 1j*eta))) for H in hamiltonian_op.GetMatrix()]) for E in get_energy_grid() ])*2
dos_theo = np.array([ 1/np.sqrt(4- E*E ) for E in get_energy_grid() ])*2*10000


my_invcos = [np.imag(H[0,0]) for H in green_func_op.GetMatrix(0.0)]
dos = np.imag(my_trace.compute( green_func_op))

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("tkagg")  # Use TkAgg backend
plt.ion()  # Turn on interactive mode

#plt.plot([x[0] for x in get_kmesh()], my_cos )
#plt.plot([x[0] for x in get_kmesh()], my_invcos_exact )
#plt.plot([x[0] for x in get_kmesh()], my_invcos, label="3" )

plt.plot(get_energy_grid(), dos_exact)
plt.plot(get_energy_grid(), dos)
plt.plot(get_energy_grid(), dos_theo)

plt.show()
input("Press Enter to exit...")  # Keeps the window open



#my_trace.compute(green_func_op)

#.compute(block_matrix=Hamitlonian.matrix)



#dos  = Trace ( GreenFunction, kgrid ).compute()

#print("Filling the hamiltonians")
#ham_k = kdens.operator( hamiltonian ); 

#print("computing eigenvalues")
#U_k   = kdens.eigenU(ham_k)
#H_new = kdens.change_basis(U_k, ham_k);

#print("diagonal elements")

#eigvals= np.array([ np.diag(h) for h in H_new]).flatten();


#Chebyshev polynomials
#emin = np.min(eigvals);
#emax = np.max(eigvals);

#alpha = 0.9;
#W = (emax+emin)/2;
#DE = emax - emin;
#eigvals = 2*alpha*(eigvals - W )/DE; 

#DOS
#M=3000;
#mu = [];

#print("KPM")

#for m in range(M):
#    T_m = np.cos( m *np.arccos(eigvals) );
#    mu.append( np.sum(T_m) );

#print(mu)
