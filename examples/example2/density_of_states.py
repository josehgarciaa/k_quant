#Run cell
#%%
import numpy as np
import matplotlib.pyplot as plt
from k_quant.system import System 

import kpm 
kpm.safe_CUTOFF= 0.99;

wann_syst = System( dimensions = (10000,1,1), w90_inp="linear_chain")
dens = kpm.Density(wann_syst, bounds=(-1,1))
dens.ComputeMoments( broadening =0.05)

plt.plot(*dens.spectral_average() )
# %%
