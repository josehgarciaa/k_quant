from numpy import sqrt, exp, dot, conj
import matplotlib.pyplot as plt
import k_quant as k
import matplotlib
matplotlib.use('TKAgg')

"""
Let us first define a Hamiltonian in momentum space
and for this purpose we choose as prototype a model for p_z electrons
 in graphene within the nearest neighbor approximation with lattice vectors
 defined in the lat_vec variable
"""
lat_vec = [ [ 1/2, sqrt(3)/2,0 ], [ 1/2,-sqrt(3)/2,0 ], [ 0,0,1 ] ]
def hamiltonian(k):
    a_0, a_1, a2  = lat_vec;
    hop = 2.8;
    f_k = hop*( 1 + exp( -1j*dot(k,a_0)) + exp( -1j*dot(k,a_1)) );
    return [ [ 0        , f_k],
             [ conj(f_k),  0 ]
            ];

#The band class requires lattice vectors and the hamiltonian function
graphene = k.bandstructure(lat_vec, hamiltonian );

#To plot a desire band-path you pass it to the class as a list of tuples
npts= 100;
bandpath = [ ("K", (1/3,2/3,0), 111 ), ("G", (0,0,0), 35) , ("M",(1/2,1/2,0),55), ("K'",(2/3,1/3,0),1) ];
graphene.set_bandpath(bandpath);

#The computing the band structure is as simmple as
bandstructure = graphene.compute_bands();
xaxis = graphene.Xaxis();
xlabels= graphene.XLabels();
for band in bandstructure:
    plt.plot(xaxis, band);

plt.gca().set_xticks(xlabels[0])
plt.gca().set_xticklabels(xlabels[1])
plt.savefig('graphene_band_structure.pdf');

sigma_x,sigma_y = [ [[0,1],[1,0]], [[0,-1j],[1j,0]] ];
bandstructure = graphene.compute_bands( proj_ops=[ sigma_x,sigma_y] );

print(bandstructure.shape)

for proj_band in bandstructure:
    band,sigma_x,sigma_y = proj_band ;
    plt.plot(xaxis,band);

plt.gca().set_xticks(xlabels[0])
plt.gca().set_xticklabels(xlabels[1])
plt.savefig('graphene_band_structure2.pdf');
