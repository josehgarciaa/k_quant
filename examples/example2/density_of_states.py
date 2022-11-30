import numpy as np
import k_quant as k

from kpm import Density

wann_syst = k.System( dimensions = (3,3,1), w90_inp="linear_chain" )

dos = Density(wann_syst, minimum_broadening = 3 );

dos.compute_moments() 


