import numpy as np

from k_quant.system import System 
from kpm import Density

wann_syst = System( dimensions = (100,1,1), w90_inp="linear_chain")
dens = Density(wann_syst)

dens.ComputeMoments( broadening =0.1)

