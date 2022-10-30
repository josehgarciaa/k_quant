A simple band structure calculation
====================================

For this example, let us first define a Hamiltonian in momentum space which we will
use to compute the band structure. In this case, we choose as prototype a model for :math:`p_z` electrons
in graphene within the nearest neighbor approximation with lattice vectors
defined in the lat_vec variable

.. code-block:: python

    from numpy import sqrt, exp, dot, conj, min,max, abs
    lat_vec = [ [ 1/2, sqrt(3)/2,0 ], [ 1/2,-sqrt(3)/2,0 ], [ 0,0,1 ] ]
    def hamiltonian(k):
        a_0, a_1, a2  = lat_vec;
        hop = 2.8;
        f_k = hop*( 1 + exp( -1j*dot(k,a_0)) + exp( -1j*dot(k,a_1)) );
        return [ [ 0        , f_k],
                [ conj(f_k),  0 ]
                ];


Then, we need to define the bandpath using as a list of tuples, consisting in 
the kpoint label, its position in the normalized Brilluoin zone, and the numbers
of point between that point and the next bandpath point. Once define, this variable
is pass to the bandstructure class.

.. code-block:: python

    #The band class requires lattice vectors and the hamiltonian function
    graphene = k.bandstructure(lat_vec, hamiltonian );

    #To plot a desire band-path you pass it to the class as a list of tuples
    npts= 100;
    bandpath = [ ("K", (1/3,2/3,0), 111 ), ("G", (0,0,0), 35) , ("M",(1/2,1/2,0),55), ("K'",(2/3,1/3,0),1) ];
    graphene.set_bandpath(bandpath);

    #The computing the band structure is as simmple as
    sigma_x,sigma_y = [ [[0,1],[1,0]], [[0,-1j],[1j,0]] ];
    bandstructure = graphene.compute_bands( proj_ops=[ sigma_x,sigma_y] );


Finally we plot the result

.. code-block:: python

    fig= plt.gcf();
    ax = plt.gca();
    xaxis = graphene.Xaxis();
    for proj_band in bandstructure:
        band,sigma_x,sigma_y = proj_band;
        z   = sigma_x;
        s   =  20*abs(z);
        plt.plot(xaxis,band, c="k");
        im  = ax.scatter(xaxis,band,s=s, c=z,cmap="coolwarm",vmin=-1, vmax=1);
    fig.colorbar(im, ax=ax);

    xlabels= graphene.XLabels();
    ax.set_xticks(xlabels[0])
    ax.set_xticklabels(xlabels[1])
    plt.savefig('proj_band_sigma_x.pdf');
    plt.show();
