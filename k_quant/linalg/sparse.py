import numpy as np
import pandas as pd
import scipy as scp
from scipy.sparse import bmat, csr_matrix, block_diag,eye,identity, diags

block_mat = bmat;
CSR       = csr_matrix;
BlockDiag = block_diag
Diag      = diags
eye       = eye
identity  = identity;



def write_csr( A, ofname ):
    dims = A.shape;
    assert dims[0]==dims[1], f"The generated operator is not a square matrix. Its shape is : {dims}";   
    rank = str(dims[0]);
    nnz  = str(A.nnz);

    
    complex2str = np.vectorize(lambda x: str(x.real)+" "+str(x.imag)); 
    integer2str = np.vectorize(lambda x: str(x));
    
    with open(ofname, "w") as f:
        data    = pd.Series(complex2str(A.data)).str.cat(sep=' ')
        indices = pd.Series(integer2str(A.indices)).str.cat(sep=' ') 
        indptr  = pd.Series(integer2str( A.indptr)).str.cat(sep=' ')
        f.write(rank+ " "+nnz+"\n"+data+"\n"+indices+"\n"+indptr+"\n")
    return ;
