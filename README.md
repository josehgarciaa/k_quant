# What is k-quant
Quantum transport calculations in k-space

https://josehgarciaa.github.io/k-quant/


## Installation/Usage:

To install the package download the code from its [repository](https://github.com/josehgarciaa/k-quant), 
then from the k-quant directoy proceed with 

```
pip install . 
```

## An example of k-quant

```python

from numpy import sqrt, exp, dot, conj
from kquant import bands

#Let us first consider the graphene Hamiltonian in momentum space
def hamiltonian(k):
    a_0 = [ 1/2, sqrt(3)/2 ];
    a_1 = [ 1/2,-sqrt(3)/2 ];
    hop = 2.8;
    f_k = hop*( 1 + exp( -1j*dot(k,a_0)) + exp( -1j*dot(k,a_0)) );
    return [ [ 0        , f_k],
             [ conj(f_k),  0 ]
            ];


```