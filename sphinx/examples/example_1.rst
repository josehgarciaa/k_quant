Calculation of the density of states
====================================

Any calculation in KQUANT uses the as basis the eigenvector of the Hamiltonian as main basis even though
the model and all operators are defined in the momentum space. In this example we will see how we harmonized
this seemly contraction. 


For this example, we use as a model a hamiltonian extracted from a wannierization. 

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

Such a model is defined in the bloch space in the base of the reciporcal and fractional coordinates
for momentum and real-space vectors respectively. Therefore, it does not contain information of the
underlying structure. The :py:class:`Lattice` class is in charge of handling structural information. 

For instance, for this particular case, we have

.. code-block:: python


    lat = Lattice( w90_in = "");

    print( lat.primitivec_vectors )
    print( lat.orbitals )
    print( lat.unit_cell )


Explain the output


By definition, a crystal lattice represents an infinite colection of point. However, KQUANT deals with systems of finite size which are not necessesarly 
periodic (see :ref:`my target` for discussion of real-space, momentum, and eigen). Therefore, we define our current system thogh the class :py:class:`System` 
which should be initialize with three dimensions

.. code-block:: python


    syst = System( dimensions=(100,100,1), w90_in = "", model="_he", deploy=True);
    
    #The lattice can be accessed through
    syst.lattice;

    #The model is accesible through
    syst.model( k= [0,0,0] );

The deploy flag is used to indicate that you want to compute all :ref:`necessary quantities`
to perform subsequent calculations. Therefore, for large systems, this may take some time, but it is performed once.

In k_quant, we don't want to commit to a particular way of computing spectral quantities. Therefore, in the module :py:module:`spectral_solvers` we
will offer different flavors for performing calculations. At this point, the most tested solvers are those based on the kernel polynomial methods
(see some relevant examples here :ref:`kpm_examples`). In this example we will use this particular solvers


.. code-block:: python

    from k_quant.spectral_solvers import kpm

    syst = System( dimensions=(100,100,1), w90_in = "", model="_he", deploy=True);

    energies = np.linspace(-1,1,100);
    dos = kpm.Density(System, broadening = 10, Op=None, energies=energies );

Then the final version of the example will


.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    import k_quant as kq
    from k_quant.spectral_solvers import kpm

    syst= kq.System( dimensions=(100,100,1), w90_in = "", model="_he", deploy=True);
    Es  = np.linspace(-1,1,100);
    dos = kpm.Density(System, broadening = 10, Op=None, energies=Es );

    plot.plot(Es,dos);









    from numpy import sqrt, exp, dot, conj, min,max, abs
    lat_vec = [ [ 1/2, sqrt(3)/2,0 ], [ 1/2,-sqrt(3)/2,0 ], [ 0,0,1 ] ]
    def hamiltonian(k):
        a_0, a_1, a2  = lat_vec;
        hop = 2.8;
        f_k = hop*( 1 + exp( -1j*dot(k,a_0)) + exp( -1j*dot(k,a_1)) );
        return [ [ 0        , f_k],
                [ conj(f_k),  0 ]
                ];



This information will be loaded 




Then, we need to define the bandpath using as a list of tuples, consisting in 
the kpoint label, its position in the normalized Brilluoin zone, and the numbers
of point between that point and the next bandpath point. Once define, this variable
is pass to the bandstructure class.

.. code-block:: python

    import k_quant as k

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
    
    import matplotlib.pyplot as plt

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
