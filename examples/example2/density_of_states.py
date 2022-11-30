import numpy as np

from system import System 


wann_syst = System( dimensions = (3,3,1), w90_inp="PtSe2" )

print("Computing eigenvalues")
Es = wann_syst.BlochEigenvalues();

print("Computing hamiltonian in eigenvalue basis")
MEs=  wann_syst.BlochToEigen( wann_syst.Hamiltonian() )

print( np.sum([ np.sum( np.abs(E1 - np.diag(ME2)) ) for E1, ME2 in zip(Es, MEs )] ) ) 



    def compute_moments(self,broadening):
        return self;
        
    def spectral_average(self,energies):
        return self;
