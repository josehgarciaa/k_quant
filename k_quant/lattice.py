from numpy.linalg import inv
from numpy import dot, pi


class Lattice():

    """ Storage and handles all structural information of a given system. 

    Args:
        primitive_vectors (array): An array containing the three primitive_vectors that spams the lattice
        orbital_positions (array): an array containing the positions and labels of all hoppings in the system

    Attributes:
        primitive_vectors (array): An array containing the three primitive_vectors that spams the lattice
        orbital_positions (array): an array containing the positions and labels of all hoppings in the system
    """
   
   
    def __init__(self, primitive_vectors= None, orbital_positions= None):
        self.primitive_vectors = primitive_vectors;
        self.orbital_positions = orbital_positions;
        pass
    
 
    def OrbitalNumber(self):
        """Returns the number of orbitals in the model

        :return: Total number of orbitals in the model
        :rtype: int
        """
        return len(self.orbital_positions);
    
    def rec2cart(self, kpoints):
        """Transforms a set of k-points in the reciprocal units in cartesian in units 
        of inverse of the primitive lattice length.

        :param kpoints: Either a k-point or list of k-points in reciprocal units 
        :type kpoints: A tuple (k_0,k_1,k_2) or a list of tuples  [(k_0^0,k_1^0,k_2^0),....,(k_0^n,k_1^n,k_2^n) ]
        :return: Either a k-point or list of k-points in cartesian in units inverse of primitive lattice length 
        :rtype: A tuple (k_0,k_1,k_2) or a list of tuples  [(k_0^0,k_1^0,k_2^0),....,(k_0^n,k_1^n,k_2^n) ]

        """
        Urc = 2*pi*inv(self.primitive_vectors).T;
        return dot( kpoints, Urc )


    