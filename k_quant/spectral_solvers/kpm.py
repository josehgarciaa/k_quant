


class Density:

    moments = None;
    Ham     = None;
    Op      = None;

    def __init__(self,system, minimum_broadening, broadening_type = "jackson",  Op=None):
        self.Ham = system.Hamiltonian();

    def compute_moments(self,broadening):
        return self;
        
    def spectral_average(self,energies):
        return self;
    
    