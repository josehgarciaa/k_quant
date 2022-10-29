
from os import wait3
import numpy as np
import numpy.typing as npt
from numpy.linalg import eigh, eigvalsh, norm
from numpy import exp, cos, sin , pi, kron

class bandstructure:
    """This class will compute the bandstructure for different inputs of models

    The __init__ method may be documented in either the class level
    docstring, or as a docstring on the __init__ method itself.

    Either form is acceptable, but the two should not be mixed. Choose one
    convention to document the __init__ method and be consistent with it.

    Note:
        Do not include the `self` parameter in the ``Args`` section.

    Args:
        msg (str): Human readable string describing the exception.
        code (:obj:`int`, optional): Error code.

    Attributes:
        msg (str): Human readable string describing the exception.
        code (int): Exception error code.

    """
       
    def __init__(self, lat_vec, H_k ):
              
        self.lat_vec   = lat_vec;
        self.bandpath  = None;
        self.H_k       = H_k;

    def set_hamiltonian_k(self, H_k):
        """_summary_

        Args:
            H_k (_type_): _description_
        """
        self.H_k = H_k;

    def hamiltonian_k(self, k):
        """_summary_

        Args:
            k (_type_): _description_

        Returns:
            _type_: _description_
        """
        return self.H_k(k);

    def set_bandpath(self,bandpath):
        """_summary_

        Args:
            bandpath (_type_): _description_
            absolute_coords (bool, optional): _description_. Defaults to False.
        """
        self.bandpath = bandpath;

    def bandpath(self):
        """_summary_

        Args:
            bandpath (_type_): _description_
            absolute_coords (bool, optional): _description_. Defaults to False.
        """
        return self.bandpath;


    def Momentum_Rec2AbsMatrix(self ):
        return 2*np.pi* np.linalg.inv(self.lat_vec);
    
    #If required, rescale to absolute value
    def toAbsoluteCoords(self,x):
             return np.dot( x, np.transpose(self.Momentum_Rec2AbsMatrix() ) );
                

    def band_kpoints(self , absolute_coords = True): 
        """Compute the k-points defined in the bandpath.  

        Args:
            absolute_coords (bool, optional): Determine whereas the k-points are given in cartesian (in units of 2pi/a) or reciprocal units. Defaults to False.

        Returns:
            ndarray: An array of k-points
        """      
        labels, segments, npoints = list(zip(*(self.bandpath)));
        segments = np.array(segments);
        last_pt  = segments[-1]; 
        segments = np.transpose( [ segments[:-1],segments[1:] ], axes=(1,0,2) )

        kpoints  = [ np.linspace( *seg, n,endpoint=False ) for seg,n in zip(segments,npoints) ] ;
        kpoints  = np.concatenate( [*kpoints,[last_pt]] , axis=0)
        #If required, rescale to absolute value
        if absolute_coords  is True :
            kpoints = np.dot( kpoints, np.transpose(self.Momentum_Rec2AbsMatrix() ) );

        return np.array(kpoints);

    def XLabels(self ):
        """Returns the x-axis labels associated with the bandpath

        Returns:
            list: _description_
        """
        Xaxis = self.Xaxis()    
        labels, segments, npoints = list(zip(*(self.bandpath)));
        xpos = np.concatenate([[0], np.cumsum( npoints)])[:-1]

        return [list(Xaxis[xpos]),labels]                
       
    def Xaxis(self ):
        """Returns a unique x-axis parameter that respect the distance between k-points

        Returns:
            ndarray: An array of k-point
        """
        Xaxis=np.cumsum(np.linalg.norm(np.diff( self.band_kpoints() ,axis=0, prepend = 0),axis=1));
        Xaxis-=Xaxis[0];#Remove the initial value.
       
        return Xaxis;
    

    def operator_k(self, operator, kpoint ):
        """_summary_

        Args:
            proj_ops (_type_): _description_

        Returns:
            _type_: _description_
        """

        #In numpy eigh, the put is a matrix where column v[:, i] 
        # is the normalized eigenvector of w[i] eigenvalue
        w, v = np.linalg.eigh( self.H_k(kpoint) );

        operator = np.array(operator(kpoint) if callable(operator) else operator, dtype=complex);
        try:
            ops_k = np.sum(np.conj(v)* operator.dot(v), axis=0) ;
            return np.real(ops_k);
        except:
            print("problem with operators_k");
                
        return None;

    
    def compute_bands(self, fermi_energy = 0.0, proj_ops = None ):
        """_summary_

        Args:
            fermi_energy (float, optional): _description_. Defaults to 0.0.
            proj_op (_type_, optional): _description_. Defaults to None.
            ax (_type_, optional): _description_. Defaults to None.
            plot_proj (bool, optional): _description_. Defaults to False.
            proj_range (_type_, optional): _description_. Defaults to None.

        Returns:
            _type_: _description_
        """

        #Compute the kpoints based on the path
        if proj_ops is None:
            bands = np.transpose( [ self.operator_k( self.H_k(kp), kp)   for kp in  self.band_kpoints() ]) - fermi_energy;
            return bands;

        ops = [self.H_k, *proj_ops];
        proj_bands = np.array( [ [ self.operator_k(op, kp) for op in  ops] for kp in  self.band_kpoints()] );
        
        return np.transpose(proj_bands, axes=(2,1,0));

def energies_inwindow(bands, energy_window):
    Emin,Emax = energy_window;
    return np.any(bands>Emin,axis=1)*np.any(bands<Emax,axis=1);  

