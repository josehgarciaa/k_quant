#Run cell
#%%

import numpy as np
import matplotlib.pyplot as plt
from k_quant.system import System 
from kpm import Density

wann_syst = System( dimensions = (4,1,1), w90_inp="linear_chain")
dens = Density(wann_syst, bounds=(-4,4))

dens.ComputeMoments( broadening =0.01)

plt.plot(*dens.spectral_average() )
# %%
