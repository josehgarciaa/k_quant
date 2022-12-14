# What is k-quant
Quantum transport calculations in k-space

https://josehgarciaa.github.io/k_quant/


![alt text](https://github.com/josehgarciaa/k_quant/blob/main/sphinx/logos/k_quant_logo.png)




## Installation/Usage:

To install the package download the code from its [repository](https://github.com/josehgarciaa/k-quant), 
then from the k-quant directoy proceed with 

```
pip install . 
```

## An example of k-quant

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

#The band class requires lattice vectors and the hamiltonian function
graphene = k.bands(lat_vec, hamiltonian );

#To plot a desire band-path you pass it to the class as a list of tuples
npts= 100;
bandpath = [ ("K", (1/3,2/3,0), npts ), ("G", (0,0,0), npts) , ("M",1/2,1/2,0) ];
graphene.set_bandpath( bandpath);

#The computing the band structure is as simmple as
bandstructure = graphene.compute_bands();
for band in banstructure():
    plt.plot(band);

#In this particular scenario, the distance between the high-symmetry point in the band-structure
#do not corresponds to the real distance within the Brilluoin Zone. This is why we provide 
#function to get the proper xaxis in the bandstructure plot
xaxis =  graphene.Xaxis();
for band in banstructure():
    plt.plot(xaxis,band);

#In several situations, it is desireable to compute the projection of orbital or spins in these bands. 
#For instance, in graphene, one would like to compute the so called pseudo-spin projected bandstructures. 
#For that case we need to define the projection operator
sigma_z = [ [1, 0],
            [0,-1]]

# then we can proceed with the calculation of the projection
#The computing the band structure is as simmple as
projections = graphene.compute_bands( projection_operators = sigma_z);
for band in banstructure():
    plt.plot(band);


```
