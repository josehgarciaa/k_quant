from k_quant.utils import create_mesh

class Grid:
    """
    Represents the geometry of a 2D mesh grid based on k-point windows.

    Attributes:
        lat_vec (np.ndarray): Lattice vectors of the system.
        kgrid (np.ndarray): The 2D mesh grid in reciprocal space.
    """


    def __init__(self, lat_vec, kdims):
        """
        Initializes the Grid object with lattice vectors and momentum grid dimensions.

        Args:
            lat_vec (np.ndarray): Lattice vectors of the system.
            kdims (tuple): A three-dimensional tuple defining the dimensions of the momentum grid.
        """
        self.lat_vec = lat_vec
        self.kdims = kdims


    def resize(self, kdims) :
        """
        Generates a 2D mesh grid based on the supercell dimension (kdims).

        Args:
            kdims (tuple): A three-dimensional tuple defining the dimensions of the momentum grid.


        """
        self.kdims = kdims
        return self
    

class Density:

    #This function creates a 2D mesh grid based on a kpoint windows kwindow=[kmin,kmax]
    # kmin is the position of the origin of the windows while kmax is the position of the
    # end of the rectangle
    #    ________ kmax
    #   |        |
    #   |        |
    #   |        |
    #   |________|
    # kmin

    def __init__(self, lat_vec, H_k ):
              
        self.lat_vec   = lat_vec;
        self.bandpath  = None;
        self.kgrid     = self.bz_grid( scdim=(1,1,1) );

    def set_scdim(self, scdim):
        return self.bz_grid(scdim=scdim );

    def transform_rec2cart(self, kpoints):
        rec2cart = 2*np.pi* np.linalg.inv(self.lat_vec).T;
        return np.dot( kpoints, rec2cart )

    def bz_grid(self, scdim=(1,1,1) ):
        meshkgrid  = np.meshgrid( *[ np.linspace(0,1,d, endpoint=False)  for d in scdim ], indexing='ij');  
        self.kgrid = np.transpose( [ x.flatten() for x in meshkgrid] );      
        return self;

    def operator(self, Op ):
        return np.apply_along_axis( Op, arr=self.kgrid, axis=1 );

    def energies_inwindow(self, bands, energy_window):
        Emin,Emax = energy_window;
        return np.any(bands>Emin,axis=1)*np.any(bands<Emax,axis=1);

    def CSR_representation( self ):
        Aop= bdiag_mat( self.Dorb, format="csr");
        Aop.eliminate_zeros();
        return Aop