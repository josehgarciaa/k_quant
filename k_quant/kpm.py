
"""
main.py
====================================
The core module of my example project
"""

import numpy as np
import k_quant.linalg.sparse  as sp
import k_quant as k

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
    """The spectral density of an operator :math:`X`.


        The spectral density is defined as :math:`\langleX\rangle(E) = {\rm Tr}[ X \delta(H-E) ]` and 
        quantifies how much of the operator X will be measured at a given en energy 
        
        In this module  :math:`\delta(H-E)` is computed using the kernel polynomial method.  



        If the class has public attributes, they may be documented here
        in an ``Attributes`` section and follow the same formatting as a
        function's ``Args`` section. Alternatively, attributes may be documented
        inline with the attribute's declaration (see __init__ method below).

        Properties created with the ``@property`` decorator should be documented
        in the property's getter method.

        Attributes
        ----------
        attr1 : str
            Description of `attr1`.
        attr2 : :obj:`int`, optional
            Description of `attr2`.

    
        A spectral density, is a quantity defined as <X> = Tr[ X delta(H-E) ]. From its definition, 
        this quantity computes how much of the operator X will be measured at a given en energy. 
        
        In the kpm module, the delta function is computed by expanding it in terms of Chebyshev polynomials
        and regularized using a kernel function.  

        
        Exceptions are documented in the same way as classes.

        The __init__ method may be documented in either the class level
        docstring, or as a docstring on the __init__ method itself.

        Either form is acceptable, but the two should not be mixed. Choose one
        convention to document the __init__ method and be consistent with it.

        Note
        ----
        Do not include the `self` parameter in the ``Parameters`` section.

        Parameters
        ----------
        msg : str
            Human readable string describing the exception.
        code : :obj:`int`, optional
            Numeric error code.

        Attributes
        ----------
        msg : str
            Human readable string describing the exception.
        code : int
            Numeric error code.
            
            
        .. _references
        References:
        [1] https://doi.org/10.1016/j.physrep.2020.12.001

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
        """Example of docstring on the __init__ method.

        The __init__ method may be documented in either the class level
        docstring, or as a docstring on the __init__ method itself.

        Either form is acceptable, but the two should not be mixed. Choose one
        convention to document the __init__ method and be consistent with it.

        Note
        ----
        Do not include the `self` parameter in the ``Parameters`` section.

        Parameters
        ----------
        param1 : str
            Description of `param1`.
        param2 : :obj:`list` of :obj:`str`
            Description of `param2`. Multiple
            lines are supported.
        param3 : :obj:`int`, optional
            Description of `param3`.

        Construct an instance of the Density class 

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
        """Generate a random phase state

        :return: _description_
        :rtype: _type_
        """
        
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
    
    #@measure_time
    def ComputeMoments(self,broadening = None):

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
        
    def spectral_average(self,energies = None):

        if energies is None:        
            energies = np.linspace(-safe_CUTOFF,safe_CUTOFF, 1000);
            densities= [np.sum([ 2*mu*np.cos(m*np.arccos(x)) for m,mu in enumerate(self.moments)])/np.sqrt(1-x**2) for x in energies];        
        return (energies, densities)
        
    
    
    def OriginalHam(self):
        print("The function OriginalHa, is not implemented yet")
        return self.Ham;