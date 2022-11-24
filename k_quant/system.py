from numpy import meshgrid, transpose,linspace, apply_along_axis
from .models.wannier import WannierSystem  
from .lattice import Lattice
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
    Ham_fun   = None;
    reciprocal_lattice = None; 

    def __init__(self, dimensions, **kwargs):

        self.dimensions= dimensions;
        self.bz_points = self.bz_points();

        #When the user submit a Wannier90 input, use it to initialice the system
        if "w90_inp" in kwargs:
            w90_fname = kwargs["90_inp"];
            syst = WannierSystem(w90_fname);
            self.lattice = Lattice( primitive_vectors = syst.primitive_vectors,
                                    orbital_positions = syst.primitive_vectors );
            self.Ham_fun = syst.hamiltonian;

        if "model" in kwargs:
            w90_fname = kwargs["model"];
            syst = WannierSystem(w90_fname);
            self.lattice = Lattice( primitive_vectors = syst.primitive_vectors,
                                    orbital_positions = syst.primitive_vectors );
            self.Ham_fun = syst.hamiltonian;


            
    def reciprocal_lattice(self):
        """Returns the points in the first Brilluoin zone that defines the reciprocal lattice.
            
            Note:
            This functions compute these points on a first call. Therefore, that first call could take significantly more time than subsequent calls.    
        
        """
        if self.reciprocal_lattice is None:
            
            dim = self.dimensions;
            meshkgrid    = meshgrid( *[ linspace(0,1,d, endpoint=False)  for d in dim ], indexing='ij');  
            self.reciprocal_lattice = transpose( [ x.flatten() for x in meshkgrid] );      

        return self.reciprocal_lattice;   
            

    def H_k(kpoint):
        """ The k-dependent hamiltonian of the system evaluted in a particular kpoint

        :param kpoint: An arbitrary kpoint within the Brilluoin zone
        :type kpoint: a tuple consisting of three real numbers
        :return: A matrix of the hamiltonian evaluted at that particular kpoint
        :rtype: Array
        """
        return self.Ham_fun(kpoint);


    def Hamiltonian(self):
        """The hamiltonian matrix in the Block diagonal sparse format the first Brilluoin zone that defines the reciprocal lattice.

        :return: A complex sparse matrix in the block-diagonal form
        """
        return apply_along_axis( self.H_k, arr=self.reciprocal_lattice, axis=1 );
            
            
    
    
