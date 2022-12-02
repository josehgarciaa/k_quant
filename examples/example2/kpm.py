import numpy as np
import k_quant.linalg.sparse  as sp
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
        self.num_kpts = system.KpointNumber();
        self.num_orbs = system.OrbitalNumber();
        
        if bounds is None:
            self.StocasticBounds();
        else:
            self.bounds = bounds;
        print("Initializing k_quant.kpm.Density will rescale the Hamiltonian spectrum : (",-safe_CUTOFF,safe_CUTOFF,")")
        self.Ham.Rescale(self.ScaleFactor(),-self.ShiftFactor());

        self.StocasticBounds();

        print("If you plan to use the original Hamiltonian call k_quant.kpm.Density.OriginalHam()")


    
    def StochasticStates(self):
        nkpt= self.num_kpts;
        neig= self.num_orbs;
        dim    = self.num_kpts * self.num_orbs;
        states = np.zeros((dim, neig), dtype=complex);       
        for n in range(neig):
            expPhi = np.exp(2j*np.pi*np.random.rand(nkpt))/np.sqrt(dim)
            states[ n::neig , n] = expPhi; #Add an exponential every neig jumps
            states[ : , n] = (self.Ham.U)@states[:,n];
        return states;

    def StocasticBounds(self):
        Rstates = self.StochasticStates();
        Lstates = np.conj(Rstates.T);

        #R2 = np.array( [ self.Ham.dot(x) for x in Rstates.T]).T;
        Rstates = self.Ham@Rstates
        
        Ham_moms = [];
        #Ham_moms2 = [];
        for n in range(4):
            Ham_moms.append( np.trace(Lstates@Rstates) );
            #Ham_moms2.append( np.sum( [np.dot(xL, xR) for xL,xR in zip(Lstates,R2.T)]) ); 
            Rstates = self.Ham@Rstates
            #R2 = np.array( [ self.Ham.dot(x) for x in R2.T]).T;
     
#        #Ham_mom are the statistical moments of the hamiltonian, ie = <H>, <H^2>, <H^3>, etc 
        Em  = np.real(Ham_moms[0]);

#        # <(H-Em)**2> = H^2 -2*Em<H> + Em**2 = H**2 - Em**2
        E2m = Ham_moms[1] - Em**2;  
        
        # <(H-Em)**4> = <H^4> - 4 <H^3>Em + 6 <H^2> Em*2 Em**2 - 4*<H>Em**3 + Em**4 
        E4m = Ham_moms[3] - 4*Ham_moms[2]*Em + \
              6*Ham_moms[1]*Em**2 - 3*Em**4; 

        W = 2.5*np.sqrt( np.real(E4m)/np.real(E2m) )
        self.bounds = ( (2*Em - W)/2 , (2*Em + W)/2 );
        
        print("Bounds computing using the stochastic approach are",self.bounds)

        return self;
    
    def ShiftFactor(self):
        return 2*safe_CUTOFF/( self.bounds[1] - self.bounds[0] )*( self.bounds[1] + self.bounds[0] );
    

    def ScaleFactor(self):
        return 2*safe_CUTOFF/( self.bounds[1] - self.bounds[0] ) ;

              
    def BroadeningToMoments( self, broadening):
        print(self.bounds)
        Emax, Emin = self.bounds;
        print("kpm.py: broadeningtomoments should be checked")
        return  int( np.pi*self.ScaleFactor()/ broadening )
    
    def ComputeMoments(self,broadening = None):

        num_mom = self.BroadeningToMoments(broadening) ;
        print("kpm.py: computing moments using broadening ",broadening, " and ", self.BroadeningToMoments(broadening) )
        moments  = np.zeros(num_mom, dtype=complex)
       
        Phi0 = self.StochasticStates();
        PhiL = np.conj(Phi0.T); 
        moments[0] = np.trace(PhiL@Phi0);     
        
        Phi1 = self.Ham@(Phi0);
        moments[1] = np.trace(PhiL@Phi1);

        for i in np.arange(2, len(moments)):
            Phi0      = 2.0 * self.Ham @ Phi1 - Phi0; 
            Phi0,Phi1 = Phi1,Phi0;
            moments[i]= np.sum(PhiL@Phi0);
            print(moments[i])
            
        return self;
        
    def spectral_average(self,energies):
        return self;
    
    
    def OriginalHam(self):
        print("The function OriginalHa, is not implemented yet")
        return self.Ham;