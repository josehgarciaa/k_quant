import numpy as np
import pandas as pd
from scipy.sparse import bmat, coo_matrix, block_diag

block_mat = bmat;
coo_mat   = coo_matrix;
bdiag_mat = block_diag


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
