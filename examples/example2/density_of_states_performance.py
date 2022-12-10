#Run cell
#%%



import numpy as np
import matplotlib.pyplot as plt
from k_quant.system import System 
import time

import kpm 
kpm.safe_CUTOFF= 1.0;


import time

t0 = time.time()


start = time.time()
wann_syst = System( dimensions = (1000,1,1), w90_inp="linear_chain")
end = time.time()
bsys_time = end - start;
print("The time building the system is",bsys_time)


start = time.time()
dens = kpm.Density(wann_syst, bounds=(-1,1))
end = time.time()
bden_time = end - start;
print("The time building the density is",bden_time)


start = time.time()
dens.ComputeMoments( broadening =0.01)
end = time.time()
bmom_time = end - start;
print("The time computing the moments is ",bmom_time)

tf = time.time()
etime = tf-t0
print("distribution of times")
print("The time building the system is",100*bsys_time/etime)
print("The time building the density is",100*bden_time/etime)
print("The time computing the moments is ",100*bmom_time/etime)


import os
times = [];
for nthread in (1,2,4,6,8,16, 20):
    os.environ["OMP_NUM_THREADS"] = str(nthread) 
    os.environ["OPENBLAS_NUM_THREADS"] = str(nthread) 
    os.environ["MKL_NUM_THREADS"] = str(nthread)    
    import numpy
    start = time.time()
    dens.ComputeMoments( broadening =0.01)
    end = time.time()
    bmom_time = end - start;
    print("The time computing the moments is ",bmom_time, "for ", nthread, "threads")
    times.append([nthread,bmom_time])

print(times)
#plt.plot(*dens.spectral_average() )
# %%
