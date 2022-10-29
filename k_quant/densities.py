




#This function creates a 2D mesh grid based on a kpoint windows kwindow=[kmin,kmax]
# kmin is the position of the origin of the windows while kmax is the position of the
# end of the rectangle
#    ________ kmax
#   |        |
#   |        |
#   |        |
#   |________|
# kmin

def create_2Dkgrid( kwindow=None, npoints=(1,1) ):
    if kwindow is None:
        kwindow = [ [0.0,0.0],[1.0,1.0]]
    kwindow= np.transpose(kwindow);
    kgrid = [ (np.linspace(kmin,kmax,p))  for (kmin,kmax),p in zip(kwindow,npoints) ];  
    return np.meshgrid(*kgrid,indexing='ij') ;


def energies_inwindow(bands, energy_window):
    Emin,Emax = energy_window;
    return np.any(bands>Emin,axis=1)*np.any(bands<Emax,axis=1);
