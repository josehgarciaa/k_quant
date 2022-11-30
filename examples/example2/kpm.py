import numpy as np

import k_quant as k
from k_quant.linalg import sparse as ksp
import copy
class Density:

    moments = None;
    Ham     = None;
    Op      = None;
    dims    = None;
    bounds  = None;
    broadening = None;
    
    def __init__(self,system,  broadening_type = "jackson",  Op=None, bounds=None):
        self.Ham    = ksp.BlockDiag(system.Hamiltonian());
        self.Umatrix= system.Umatrix();
        self.dims   = system.Hamiltonian().shape;
        
        if bounds is None:
            self.StocasticBounds();
        
    
    def StochasticStates(self):
        nkpt= self.dims[0];
        neig= self.dims[1];

        op       = 'ijk,ikl->ijl';
        axis_dot = np.einsum
        states = np.apply_along_axis(np.diag,  arr=np.exp(2j*np.pi*np.random.rand(nkpt,neig)), axis=1 )/np.sqrt(nkpt*neig);
        states = axis_dot( op, self.Umatrix, states );
        states = np.transpose( states, (1,0,2) ).reshape(neig, neig*nkpt )

        return states;

    def StocasticBounds(self):
        Lstates = self.StochasticStates();
        Rstates = copy.copy(Lstates);
        
        Ham_moms = [];
        for n in range(4):
            Rstates = [ self.Ham.dot(x) for x in Rstates];
            Ham_moms.append( np.sum( [np.vdot(xL, xR) for xL,xR in zip(Lstates,Rstates)]) ); 
        
        #Ham_mom are the statistical moments of the hamiltonian, ie = <H>, <H^2>, <H^3>, etc 
        Em  = np.real(Ham_moms[0]);

        # <(H-Em)**2> = H^2 -2*Em<H> + Em**2 = H**2 - Em**2
        E2m = Ham_moms[1] - Em**2;  
        
        # <(H-Em)**4> = <H^4> - 4 <H^3>Em + 6 <H^2> Em*2 Em**2 - 4*<H>Em**3 + Em**4 
        E4m = Ham_moms[3] - 4*Ham_moms[2]*Em + \
              6*Ham_moms[1]*Em**2 - 3*Em**4; 

        W = 2.5*np.sqrt( np.real(E4m)/np.real(E2m) )
        self.bounds = ( (2*Em + W)/2 , (2*Em - W)/2 );
        return self;
              
    def BroadeningToMoments( self, broadening):
        Emax, Emin = self.bounds;
        return  int( np.pi/ ( broadening/(Emax-Emin)) )
    
    def ComputeMoments(self,broadening = None):
        num_mom = self.BroadeningToMoments(broadening) ;
        print("computing moments using broadening ",broadening, " and ", self.BroadeningToMoments(broadening) )
        Lstates = self.StochasticStates();
        Rstates = copy.copy(Lstates);
        
        moments  = np.zeros(num_mom, dtype=complex)
        
        Phi0 = copy.copy(Lstates);
        Phi1 = [ self.Ham.dot(x) for x in Rstates];
        
        moments = [];
        for mu in range(num_mom):
            Phi0 = [ self.Ham.dot(x) - Phi0 for x in zip(Phi1,Phi0) ];
            Phit  = Phi1;
            Phi1  = Phi0;
            Phi0  = Phit;
        
        
       # print([np.vdot(x, x) for x in states ] )
       # print(np.sum( [np.vdot(x, self.Ham.dot(x)) for x in states]) )
        return self;
        
    def spectral_average(self,energies):
        return self;
    
    