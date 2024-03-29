import numpy as np
import k_quant.linalg.sparse as sp

class Operator:
    """ An operator is a class that can apply a lineal transformation onto a state vector \Psi. This object is representable 
        as an hermitian matrix of dimension D defined in a given basis. 
        
        In k_quant, all calculations are performed in the eigenvector basis. However, the operator
        may be defined in a different basis provided the change-of-basis matrix U is also provided.     


    Attributes
    ----------
    name : str
        first name of the person
    surname : str
        family name of the person
    age : int
        age of the person

    Methods
    -------
    info(additional=""):
        Prints the person's name and age.        

    Returns:
        _type_: _description_
    """

    mat_rep= None;
    basis  = None;
    U      = None;
    dim    = None;
    
    def __init__(self, mat_rep=None, basis_info=None):
        """
        Constructs an arbitrary operator
        
        Parameters
        ----------
            mat_rep : object
                The matrix representation of the operator. This object should be transformable into a either a BDIAG or CSR matrices
                defined in the linalg.sparse module

            basis : dictionary
                tuple tha contain as first element the name of the basis and as second the change-of-basis transformation U. 
        """
        self.basis_info = ("eigen", None);
        if basis_info is not None:
            self.basis, self.U = basis_info
            if self.U is not None:
                try:
                    self.U = sp.BlockDiag(self.U);
                except:
                    try:
                        self.U = sp.CSR(self.U);
                    except:
                        print("The matrix representation in operators is not mapable neither into CSR nor BlockDiag objects in linalg.sparse module")
       
        if mat_rep is not None:
            try:
                self.mat_rep = sp.BlockDiag(mat_rep);               
            except:
                try:
                    self.mat_rep = sp.CSR(mat_rep);
                except:
                    print("The matrix representation in operators is not mapable neither into CSR nor BlockDiag objects in linalg.sparse module")
                  

    def set_constant(self, Dorb ):
        self.Dorb = Dorb;
        return self

    def set_supercell(self, dims ):
        assert self.dims[0]==self.dims[1], f"The generated operator is not a square matrix. Its shape is : {dims}";   
        self.dims= dims
        nc       = np.prod(self.dims);
        self.Dorb= np.broadcast_to( self.Dorb, [nc,*self.Dorb.shape]).astype(complex);   
        return self


    def CSR_representation( self ):
        Aop= sp.bdiag_mat( self.Dorb, format="csr");
        Aop.eliminate_zeros();
        return Aop
    
    def __matmul__(self, A):
        return self.mat_rep@A

    def __rmatmul__(self, A):
        return A@self.mat_rep


    def __mul__(self,a):
        return self.mat_rep*a;

    def __rmul__(self,a):
        return a*self.mat_rep;
           

    def dot(self,x):
        return self.mat_rep.dot(x);
   
    def Rescale(self, a,b):
        dim0,dim1 = self.mat_rep.shape
        one = sp.identity(dim0,  dtype=complex, format= self.mat_rep.format);
        self.mat_rep = a*self.mat_rep + b*one
        return self;
