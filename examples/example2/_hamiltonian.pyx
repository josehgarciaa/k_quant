import cython
from setuptools import Extension, setup
from Cython.Build import cythonize
ext_modules = [
    Extension(
        "hamiltonian",
        ["_hamiltonian.py"],
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-fopenmp'],
    )
]

setup(
    name='hello-parallel-world',
    ext_modules=cythonize(ext_modules, compiler_directives={'language_level' : "3"} ),
)
from cython.parallel import prange


def CHamiltonian(self,k):
        """
        Calculates the Hamilton matrix for a given k-point or list of
        k-points.

        Parameters
        ----------
        k :
            The k-point at which the Hamiltonian is evaluated. If a list
            of k-points is given, the result will be the corresponding
            list of Hamiltonians.
        """

        print( "Hamiltonians")

        print("cython")
        i = cython.declare(cython.int)
        n = cython.declare(cython.int, 30)
        sum = cython.declare(cython.int, 0)

        for i in prange(n, nogil=True):
            sum += i

        print(sum)

                
        return sum;
    
