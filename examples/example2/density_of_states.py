import numpy as np

from system import System 


wann_syst = System( dimensions = (100,1,1), w90_inp="linear_chain")

print(wann_syst.Hamiltonian() )

