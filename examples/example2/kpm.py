import numpy as np
import k_quant.linalg.sparse  as sp
import k_quant as k

safe_CUTOFF = 0.95;
class Density:
    """ A spectral density, is a quantity defined as <X> = Tr[ X delta(H-E) ]. From its definition, 
        this quantity computes how much of the operator X will be measured at a given en energy. 
        
        In the kpm module, the delta function is computed by expanding it in terms of Chebyshev polynomials
        and regularized using a kernel function.  

        
    """
    
    moments = None;
    Ham     = None;
    Op      = None;
    dims    = None;
    bounds  = None;
    broadening = None;
    num_kpts= None;
    num_orbs= None;
    
    def __init__(self,system,  kernel = "jackson",  Op=None, bounds=None):
        """Construct an instance of the Density class 

        Args:
            system (object): A system as defined in the system module. 
            kernel (str, optional): The choice of kernel for the regularization. Defaults to "jackson".
            Op (Operator, optional): An operator as described in the operator module.
            bounds (tuple, optional): A tuple defining the spectral bound.

        Note:
            Using this module, will result in rescaling the hamiltonian operator defined in system. For returning it 
            to the original one, please call the method self.OriginalHam().
        """

        self.Ham    = system.Hamiltonian();
        self.num_kpts = system.KpointNumber();
        self.num_orbs = system.OrbitalNumber();
        
        if bounds is None:
            self.StocasticBounds();
        else:
            self.bounds = bounds;
        
        self.Ham.Rescale(self.ScaleFactor(),-self.ShiftFactor());

    
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
        
        print("Stochastic Bonds not implemented. Please define a bound")

        return self;
    
    def ShiftFactor(self):
        return safe_CUTOFF/( self.bounds[1] - self.bounds[0] )*( self.bounds[1] + self.bounds[0] );
    

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
        self.moments  = np.zeros(num_mom, dtype=complex)
       
        Phi0 = self.StochasticStates();
        PhiL = np.conj(Phi0.T); 
        self.moments[0] = 0.5*np.trace(PhiL@Phi0);     
        
        Phi1 = self.Ham@(Phi0);
        self.moments[1] = np.trace(PhiL@Phi1);

        for i in np.arange(2, len(self.moments)):
            Phi0= 2.0 * self.Ham @ Phi1 - Phi0; 
            Phit=Phi1; Phi1 = Phi0; Phi0=Phit;
            self.moments[i]= np.sum(PhiL@Phi0);
            
        self.moments = np.real(self.moments)
        return self;
        
    def spectral_average(self,energies = None):

        if energies is None:        
            energies = np.linspace(-safe_CUTOFF,safe_CUTOFF, 1000);
            densities= [np.sum([ mu*np.cos(m*np.arccos(x))/np.sqrt(1-x**2) for m,mu in enumerate(self.moments)]) for x in energies];        
        return (energies, densities)
        
    
    
    def OriginalHam(self):
        print("The function OriginalHa, is not implemented yet")
        return self.Ham;