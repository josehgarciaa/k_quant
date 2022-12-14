import numpy as np

safe_CUTOFF = 0.95;

import time
#def measure_time(function):
#    def wrapper(arg1, broadening):
#        print("computing moments")
#        start = time.time()
#        function(arg1, broadening)
#        end = time.time()
#        etime = end - start;
#        print("computing moments took ",etime, "seconds")
#    return wrapper


class Density:
    r"""The spectral density of an operator :math:`X`.


        The spectral density is defined as :math:`X(E) = {\rm Tr}[ \hat{X} \delta(H-E) ]` and 
        quantifies how much of the operator X will be measured at a given en energy .
        
        In this module  :math:`\delta(H-E)` is computed using the kernel polynomial method (`KPM`_).  

        Parameters
        ----------

        syst : :class:`k_quant.system.System`
            A properly initialized system.
        bounds: :obj:`str`
            The energy bounds (in the same units as the Hamiltonian) used to rescaled the hamiltonian spectrum to the (-1,1) interval. 
        kernel : :obj:`str`, optional
            A string indicating the kernel to be used for the regularization. The options are: "jackson, lorentz". 
            For any other option will use default: Jackson
        X : :class:`k_quant.operator`, optional
            The operator used to compute the density. When none submitted will use identity
                    
        Note
        ----
            Using this module, will result in rescaling the hamiltonian operator defined in system. For returning it 
            to the original one, please call the method self.OriginalHam().       
    """
    
    
    moments = None;
    Ham     = None;
    Op      = None;
    dims    = None;
    bounds  = None;
    broadening = None;
    num_kpts= None;
    num_orbs= None;
    
    def __init__(self,syst,bounds,  kernel = "jackson",  X=None ):

        self.Ham    = syst.Hamiltonian();
        self.num_kpts = syst.KpointNumber();
        self.num_orbs = syst.OrbitalNumber();
        
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
    
    def ShiftFactor(self) -> float: 
        """  Returns the shift used in the rescaling operation
        
        """
        return safe_CUTOFF/( self.bounds[1] - self.bounds[0] )*( self.bounds[1] + self.bounds[0] );
    

    def ScaleFactor(self) -> float:
        """  Returns the scale factor used in the rescaling operation
        
        """
        return 2*safe_CUTOFF/( self.bounds[1] - self.bounds[0] ) ;

              
    def BroadeningToMoments( self, broadening) -> int:
        """  Returns the number of moments for a given broadening and kernel

        Returns
        -------
            The moments are computed following the recipes in `KPM`_
            
        """

        Emax, Emin = self.bounds;
        print("kpm.py: broadeningtomoments should be checked")
        return  int( np.pi*self.ScaleFactor()/ broadening )
    
    #@measure_time
    def ComputeMoments(self,broadening) -> list:
        """Returns the chebyshev moments required to achieve a given broadening

        Parameters
        ----------
        broadening : :obj:`float`
            The value of the broadening in the same units as the Hamiltonian

        """

        num_mom = self.BroadeningToMoments(broadening) ;
        print("kpm.py: computing moments using broadening ",broadening, " and ", self.BroadeningToMoments(broadening) )
        self.moments  = np.zeros(num_mom, dtype=complex)
       
        Phi0 = self.StochasticStates();
        PhiL = np.conj(Phi0.T); 
        self.moments[0] = 0.5*np.sum(PhiL@Phi0);     
     
        Phi1 = self.Ham@(Phi0);
        self.moments[1] = np.sum(PhiL@Phi1);

        for m in np.arange(2, len(self.moments)):
            Phi0 = 2.0 * self.Ham @ Phi1 - Phi0;
            self.moments[m]= np.sum(PhiL@Phi0);
            Phit=Phi1; Phi1 = Phi0; Phi0=Phit;
            
        self.moments = np.real(self.moments)
          
        return self;
        
    def SpectralAverage(self,energies) -> list: 
        """Returns the spectral average of :math:`X` in a set of energies

        Parameters
        ----------
        energies : :obj:`list`
            The set of energies where the average will be computed

        """
        
        if energies is None:        
            energies = np.linspace(-safe_CUTOFF,safe_CUTOFF, 1000);
            densities= [np.sum([ 2*mu*np.cos(m*np.arccos(x)) for m,mu in enumerate(self.moments)])/np.sqrt(1-x**2) for x in energies];        
        return (energies, densities)
        
    
    
    def OriginalHam(self):
        print("The function OriginalHa, is not implemented yet")
        return self.Ham;