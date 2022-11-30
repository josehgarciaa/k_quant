import numpy as np
from numpy import meshgrid, transpose,linspace, apply_along_axis
from k_quant.models.wannier import WannierSystem  
from k_quant.lattice import Lattice
from k_quant.linalg.sparse import BlockDiag


class System():
    """ Storage and handle the momentum-space model and structural information of the system  

    Args:
        dimensions (tuple): A three integer tuple that defines the dimensions of the system in terms of repeating unit cells. 
        
    :Keyword Arguments:

        The arguments are presented in order of preference, i.e, w90_inp will be used over model. 
        
        w90_inp (str): The filename of a Wannier90 input file. Calling this argument will look for the label_hr.dat and label.xyz files in the same directory

    Attributes:
        dimensions (tuple): A three integer tuple that defines the dimensions of the system in terms of repeating unit cells.

    """
    
    dimensions= None;
    lattice   = Lattice();
    ham_fun   = None;
    rec_lat   = None; 
    ham_op    = None;
    Uop       = None;
    basis     = "bloch"

    def __init__(self, dimensions, **kwargs):

        self.dimensions= dimensions ;
        self.rec_lat   = self.ReciprocalLattice();

        #When the user submit a Wannier90 input, use it to initialice the system
        if "w90_inp" in kwargs:
            w90_fname = kwargs["w90_inp"];
            syst = WannierSystem(label=w90_fname);
            self.lattice = Lattice( primitive_vectors = syst.primitive_vectors,
                                    orbital_positions = syst.primitive_vectors );
            self.ham_fun = syst.Hamiltonian;

        if "model" in kwargs:
            w90_fname = kwargs["model"];
            syst = WannierSystem(w90_fname);
            self.lattice = Lattice( primitive_vectors = syst.primitive_vectors,
                                    orbital_positions = syst.primitive_vectors );
            self.ham_fun = syst.hamiltonian;


            
    def ReciprocalLattice(self):
        """Returns the points in the first Brilluoin zone that defines the reciprocal lattice.
            
            Note:
            This functions compute these points on a first call. Therefore, that first call could take significantly more time than subsequent calls.    
        
        """
        if self.rec_lat is None:
            
            dim = self.dimensions;
            meshkgrid    = meshgrid( *[ linspace(0,1,d, endpoint=False)  for d in dim ], indexing='ij');  
            self.rec_lat = transpose( [ x.flatten() for x in meshkgrid] );      

        return self.rec_lat;   
            

    def H_k(self, kpoint):
        """ The k-dependent hamiltonian of the system evaluted in a particular kpoint

        :param kpoint: An arbitrary kpoint within the Brilluoin zone
        :type kpoint: a tuple consisting of three real numbers
        :return: A matrix of the hamiltonian evaluted at that particular kpoint
        :rtype: Array
        """
        return self.ham_fun(kpoint);


    def Umatrix(self):
        """The hamiltonian matrix in the Block diagonal sparse format the first Brilluoin zone that defines the reciprocal lattice.

        :return: A complex sparse matrix in the block-diagonal form
        """
        return apply_along_axis( self.U_k, arr=self.rec_lat, axis=1 );
            

    def Hamiltonian(self):
        """The hamiltonian matrix in the Block diagonal sparse format the first Brilluoin zone that defines the reciprocal lattice.

        :return: A complex sparse matrix in the block-diagonal form
        """
        if self.ham_op is None:
            self.ham_op = apply_along_axis( self.H_k, arr=self.rec_lat, axis=1 );

        return self.ham_op;

    def Umatrix(self):
        if self.Uop is None:
            self.Uop = np.linalg.eigh( self.Hamiltonian() )[1];
        
        return self.Uop;


    def BlochToEigen(self, blochOP):
        UR = self.Umatrix();
        UL = np.conj( np.transpose( UR, (0,2,1) ) );
        op       = 'ijk,ikl->ijl';
        axis_dot = np.einsum        
        return axis_dot(op, UL, axis_dot(op, blochOP, UR ) );
    
    def BlochEigenvalues(self):
        return np.linalg.eigvalsh( self.Hamiltonian() )
        
        
            
    
    
