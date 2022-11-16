from .base_operator import Operator
from ..sparse_matrices import  bdiag_mat
import numpy as np

class Density(Operator):
    def __init__(self):
        self.dims = (1,1,1);
        pass
    
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
        Aop= bdiag_mat( self.Dorb, format="csr");
        Aop.eliminate_zeros();
        return Aop
