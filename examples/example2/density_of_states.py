#Run cell
#%%
import numpy as np
import matplotlib.pyplot as plt
from k_quant.system import System 

<<<<<<< HEAD
wann_syst = System( dimensions = (100000,1,1), w90_inp="linear_chain")
dens = Density(wann_syst, bounds=(-4,4))

dens.ComputeMoments( broadening =0.0001)
=======
import kpm 
kpm.safe_CUTOFF= 0.99;

wann_syst = System( dimensions = (10000,1,1), w90_inp="linear_chain")
dens = kpm.Density(wann_syst, bounds=(-1,1))
dens.ComputeMoments( broadening =0.05)
>>>>>>> f65c9d1085fa6bc85a10a4f9f61f457e4e9a080e

plt.plot(*dens.spectral_average() )
# %%
