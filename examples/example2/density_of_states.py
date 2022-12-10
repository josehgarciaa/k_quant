#Run cell
#%%



import numpy as np
import matplotlib.pyplot as plt
from k_quant.system import System 
import time

import kpm 
kpm.safe_CUTOFF= 1.0;

wann_syst = System( dimensions = (1000,1,1), w90_inp="linear_chain")
dens = kpm.Density(wann_syst, bounds=(-1,1))

dens.ComputeMoments( broadening =0.1)

#plt.plot(*dens.spectral_average() )
# %%
