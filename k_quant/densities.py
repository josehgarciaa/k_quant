import numpy as np

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