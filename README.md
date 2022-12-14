

# What is 

![alt text](https://github.com/josehgarciaa/k_quant/blob/main/sphinx/logos/k_quant_logo.png)

Is a package to performe efficient quantum transport in disordered crystalline systems. 

## Capabilities 

[KQuant](https://josehgarciaa.github.io/k_quant/) currently allows to compute

- The conductivity tensor

- The momentum relaxation times

- Spin Hall conductivities

- Nonequilibrium spin-densities

- Spin-orbit torques

- Spin relaxation lengths



## Inputs

We are currently compatible with the following model-building tools:

- [Wannier90](http://www.wannier.org/)

- [PAOflow](http://www.aflowlib.org/src/paoflow/)

- [Kwant](https://kwant-project.org/)


# Installation/Usage:

To install just perform the following commands
```
git clone https://github.com/josehgarciaa/k_quant.git
cd k_quant
pip install . 
```
# The density of states using KQuant

```python

import matplotlib.pyplot as plt
from k_quant.system import System 

import kpm 
kpm.safe_CUTOFF= 0.99;

wann_syst = System( dimensions = (10000,1,1), w90_inp="linear_chain")
dens = kpm.Density(wann_syst, bounds=(-1,1))
dens.ComputeMoments( broadening =0.05)

plt.plot(*dens.spectral_average() )

```

## Projected bandstructures in KQUANT


```python
from numpy import sqrt, exp, dot, conj
import kquant as k

"""
Let us first define a Hamiltonian in momentum space
and for this purpose we choose as prototype a model for p_z electrons
 in graphene within the nearest neighbor approximation with lattice vectors
 defined in the lat_vec variable
"""
lat_vec = [ [ 1/2, sqrt(3)/2 ], [ 1/2, sqrt(3)/2 ] ]
def hamiltonian(k):
    a_0, a_1  = lat_vec;
    hop = 2.8;
    f_k = hop*( 1 + exp( -1j*dot(k,a_0)) + exp( -1j*dot(k,a_1)) );
    return [ [ 0        , f_k],
             [ conj(f_k),  0 ]
            ];

graphene = k.bands(lat_vec, hamiltonian );

npts= 100;
bandpath = [ ("K", (1/3,2/3,0), npts ), ("G", (0,0,0), npts) , ("M",1/2,1/2,0) ];
graphene.set_bandpath( bandpath);

sigma_z = [ [1, 0],
            [0,-1]]

projections = graphene.compute_bands( projection_operators = sigma_z);
for band in banstructure():
    plt.plot(band);

```