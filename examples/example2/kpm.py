


class Density:

    moments = None;
    Ham     = None;
    Op      = None;
    dims     = None;

    def __init__(self,system, minimum_broadening, broadening_type = "jackson",  Op=None):
        self.Ham    = system.Hamiltonian();
        self.Umatrix= system.Umatrix();
        self.dims   = ( system.Lattice().OrbitalNumber()  );
    def compute_moments(self,broadening = None):
        print("computing moments", self.Ham.shape )
        print(self.dims)
        return self;
        
    def spectral_average(self,energies):
        return self;
    
    