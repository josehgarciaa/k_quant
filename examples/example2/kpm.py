import numpy as np

import k_quant as k

safe_CUTOFF = 0.95;
class Density:

    moments = None;
    Ham     = None;
    Op      = None;
    dims    = None;
    bounds  = None;
    broadening = None;
    num_kpts= None;
    num_orbs= None;
    
    def __init__(self,system,  broadening_type = "jackson",  Op=None, bounds=None):
        self.Ham    = system.Hamiltonian();
        #self.Umatrix= system.Umatrix();
        self.num_kpts = system.KpointNumber();
        self.num_orbs = system.OrbitalNumber();
        
        if bounds is None:
            self.StocasticBounds();
        else:
            self.bounds = bounds;
        print("Initializing k_quant.kpm.Density will rescale the Hamiltonian spectrum : (",-safe_CUTOFF,safe_CUTOFF,")")
        print("If you plan to use the original Hamiltonian call k_quant.kpm.Density.OriginalHam()")


    
    def StochasticStates(self):
        nkpt= self.num_kpts;
        neig= self.num_orbs;

        op       = 'ijk,ikl->ijl';
        axis_dot = np.einsum
        states = np.apply_along_axis(np.diag,  arr=np.exp(2j*np.pi*np.random.rand(nkpt,neig)), axis=1 )/np.sqrt(nkpt*neig);
        
        print("kpm.py:  Stochastic state should be made compatible with operator")
        #states = axis_dot( op, self.Umatrix, states );
        states = np.transpose( states, (0,2, 1) ).reshape(neig*nkpt, neig )

        return states;

    def StocasticBounds(self):
        Rstates = self.StochasticStates();
        Lstates = np.conj(Rstates.T);

        print("kpm.py:  Stochastic state should be made compatible with operator")
        

        R1 = self.Ham@Rstates
        R2 = np.array( [ self.Ham.dot(x) for x in Rstates.T]).T;
        
        Ham_moms1 = [];
        Ham_moms2 = [];
        for n in range(4):
            Ham_moms1.append( np.trace(Lstates@R1) );
            Ham_moms2.append( np.sum( [np.dot(xL, xR) for xL,xR in zip(Lstates,R2.T)]) ); 
            R1 = self.Ham@R1
            R2 = np.array( [ self.Ham.dot(x) for x in R2.T]).T;


        print( np.array(Ham_moms1)-np.array(Ham_moms2) )

#        
#        #Ham_mom are the statistical moments of the hamiltonian, ie = <H>, <H^2>, <H^3>, etc 
#        Em  = np.real(Ham_moms[0]);
#
#        # <(H-Em)**2> = H^2 -2*Em<H> + Em**2 = H**2 - Em**2
#        E2m = Ham_moms[1] - Em**2;  
        
        # <(H-Em)**4> = <H^4> - 4 <H^3>Em + 6 <H^2> Em*2 Em**2 - 4*<H>Em**3 + Em**4 
#        E4m = Ham_moms[3] - 4*Ham_moms[2]*Em + \
#              6*Ham_moms[1]*Em**2 - 3*Em**4; 

#        W = 2.5*np.sqrt( np.real(E4m)/np.real(E2m) )
#        self.bounds = ( (2*Em + W)/2 , (2*Em - W)/2 );
        self.bounds=(-1, 1)
        print("kpm.py: Stochastic bounds should be made compatible with operator")
        return self;
    
    def BandCenter(self):
        return ( self.bounds[1] + self.bounds[0] )/2;
    

    def Scale_Factor(self):
        return 2*safe_CUTOFF/( self.bounds[1] - self.bounds[0] ) ;

              
    def BroadeningToMoments( self, broadening):
        print(self.bounds)
        Emax, Emin = self.bounds;
        print("kpm.py: broadeningtomoments should be checked")
        return  int( np.pi*self.Scale_Factor()/ broadening )
    
    def ComputeMoments(self,broadening = None):
        num_mom = self.BroadeningToMoments(broadening) ;
       
       # print("kpm.py: computing moments using broadening ",broadening, " and ", self.BroadeningToMoments(broadening) )
       
       # moments  = np.zeros(num_mom, dtype=complex)


       # Phi0 = self.StochasticStates();
       # PhiL = np.conjugate(copy.copy(Phi0).T); 

       # moments[0] = np.sum(PhiL@Phi0);
       # Phi1 = self.Ham.dot(Phi0);
       # moments[1] = np.sum(PhiL@Phi1);

       # for i in np.arange(2, len(moments)):
       #     Phi0 = self.Ham.LinearT(2.0, Phi1,-1.0, Phi0 ); 
       #     Phi0,Phi1 = Phi1,Phi0; print("Swap in kpm.py ComputeMoments should be checked")
       #     moments[i] = np.sum(PhiL@Phi0);
       #     
        return self;
        
    def spectral_average(self,energies):
        return self;
    
    
    def OriginalHam(self):
        print("The function OriginalHa, is not implemented yet")
        return self.Ham;