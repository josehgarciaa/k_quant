from numpy import sqrt, exp, dot, conj
import matplotlib.pyplot as plt
import k_quant as k

"""
Let us first define a Hamiltonian in momentum space
and for this purpose we choose as prototype a model for p_z electrons
 in graphene within the nearest neighbor approximation with lattice vectors
 defined in the lat_vec variable
"""
lat_vec = [ [ 1/2, sqrt(3)/2,0 ], [ 1/2, sqrt(3)/2,0 ] ]
def hamiltonian(k):
    a_0, a_1  = lat_vec;
    hop = 2.8;
    f_k = hop*( 1 + exp( -2j*dot(k,a_0)) + exp( -1j*dot(k,a_1)) );
    return [ [ 0        , f_k],
             [ conj(f_k),  0 ]
            ];

#The band class requires lattice vectors and the hamiltonian function
graphene = k.bandstructure(lat_vec, hamiltonian );

#To plot a desire band-path you pass it to the class as a list of tuples
npts= 100;
bandpath = [ ("K", (1/3,2/3,0), npts ), ("G", (0,0,0), npts) , ("M",(1/2,1/2,0),1), ("X",(1/2,0,0),1) , ("Y",(1/3,0,0),1) ];
graphene.set_bandpath( bandpath);

#The computing the band structure is as simmple as
print( graphene.band_kpoints() );

#bandstructure = graphene.compute_bands();
#print(bandstructure.shape)

#for band in bandstructure:
#    plt.plot(band);