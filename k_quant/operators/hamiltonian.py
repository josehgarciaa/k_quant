from .base_operator import Operator
from ..linalg import  sparse as sp
import numpy as np

class Hamiltonian(Operator):
    """The Hamiltonian is a form of operator which must be defined in the momentum space
    and has as one of its main properties the capacity to generate the U matrices of the system

    """
     
    def __init__(self, mat_rep=None):
        """
        Constructs the Hamiltonian operator.

        Parameters
        ----------
            mat_rep : object
                The matrix representation of the operator. This object should be transformable into BDIAG defined in the linalg.sparse module
        """
        
        if mat_rep is not None:      
            U = self.Umatrix(mat_rep);

        basis = ("bloch", U);
        super().__init__(mat_rep, basis = basis);


    def Umatrix(self, bloch_op):
        """Compute the change-of-basis matrix from a bloch basis matrix to the eigenvector basis

        Returns:
            _type_: _description_
        """

        return np.linalg.eigh( bloch_op)[1];


    def BlochEigenvalues(self):
        print("BlochEigenvalues Not implemented yet")
        return None;

    def BlochToEigen(self, blochOP):
        print("BlochToEigen Not implemented yer")
#        UR = self.Umatrix();
#        UL = np.conj( np.transpose( UR, (0,2,1) ) );
#        op       = 'ijk,ikl->ijl';
#        axis_dot = np.einsum        
#        return axis_dot(op, UL, axis_dot(op, blochOP, UR ) );

    
