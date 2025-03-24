"""
Module for optimizing user-defined Hamiltonian functions using Numba.
"""

import inspect
import dis
from types import FrameType
from typing import Any, Callable, Set

import numpy as np
from numba import njit, prange


def get_globals_used(func: Callable) -> Set[str]:
    """
    Analyze the given function to identify all accessed global variables.

    Parameters:
        func (Callable): The function to analyze.

    Returns:
        Set[str]: A set containing the names of global variables used.
    """
    bytecode = dis.Bytecode(func)
    global_vars = {instr.argval for instr in bytecode if instr.opname == "LOAD_GLOBAL"}
    return global_vars


def make_contiguous(obj: Any) -> Any:
    """
    Convert arrays and lists to contiguous NumPy arrays.

    Parameters:
        obj (Any): The object to convert.

    Returns:
        Any: The contiguous array if conversion was applicable, otherwise the original object.
    """
    if isinstance(obj, np.ndarray):
        return np.ascontiguousarray(obj)
    if isinstance(obj, list):
        return np.ascontiguousarray(np.array(obj))
    return obj


def optimize_k_hamiltonian(user_func: Callable) -> Callable:
    """
    JIT compile the user function with automatic detection of global variables.

    Parameters:
        user_func (Callable): The user-defined Hamiltonian function.

    Returns:
        Callable: The JIT-compiled version of the user function.
    """
    # Extract the source code of the function as defined by the user.
    source_code = inspect.getsource(user_func)
    global_vars = get_globals_used(user_func)

    # Dynamically construct the namespace for the JIT-compiled function.
    namespace = {"np": np}

    # Use the caller's globals by accessing the previous frame.
    caller_frame: FrameType = inspect.currentframe().f_back
    globals_dict = caller_frame.f_globals

    for var in global_vars:
        if var in globals_dict:
            value = globals_dict[var]
            namespace[var] = make_contiguous(value)
        else:
            print(f"Warning: Global variable '{var}' not found in the provided globals dictionary.")

    print(f"Namespace for JIT compilation: {list(namespace.keys())}")

    # Execute the source code in the isolated namespace.
    exec(source_code, namespace)

    # Retrieve the function (should be an exact copy) and JIT compile it.
    compiled_func = namespace[user_func.__name__]
    jit_compiled_func = njit(compiled_func, fastmath=True)

    # Test compilation immediately to catch errors early.
    dummy_k = np.ascontiguousarray(np.zeros(3, dtype=np.float64))
    jit_compiled_func(dummy_k)

    return jit_compiled_func


@njit(parallel=True, fastmath=True)
def numba_compute_hamiltonian(points: np.ndarray, jit_ham_func: Callable) -> np.ndarray:
    """
    Compute Hamiltonians for multiple points using a Numba-optimized Hamiltonian function.

    Parameters:
        points (np.ndarray): An array of points with shape (D, N).
        jit_ham_func (Callable): The JIT-compiled Hamiltonian function.

    Returns:
        np.ndarray: An array of Hamiltonians with shape (D, 2, 2).
    """
    D = points.shape[0]
    ham_shape = jit_ham_func(points[0]).shape
    ham = np.zeros( (D,ham_shape[0],ham_shape[1]), dtype=np.complex128)
    
    ham = np.zeros(ham_shape, dtype=np.complex128)
    for idx in prange(D):
        ham[idx] = jit_ham_func(points[idx])
    return ham


def python_compute_hamiltonian(points: np.ndarray, ham_func: Callable) -> np.ndarray:
    """
    Compute Hamiltonians for a set of points using pure Python.

    Parameters:
        points (np.ndarray): An array of points with shape (D, N).
        ham_func (Callable): The Hamiltonian function (not JIT-compiled).

    Returns:
        np.ndarray: An array of Hamiltonians with shape (D, 2, 2).
    """
    try:
        D = points.shape[0]
        ham_shape = ham_func(points[0]).shape
        ham = np.zeros( (D,ham_shape[0],ham_shape[1]), dtype=np.complex128)
        for idx in range(D):
            ham[idx] = ham_func(points[idx])
        return ham
    except Exception as e:
        print(f"Error in Hamiltonian computation (Python): {e}")
        raise


def compute_hamiltonian(points: np.ndarray, ham_func: Callable) -> np.ndarray:
    """
    Automatically choose between the Numba-optimized and pure Python implementation
    of Hamiltonian computation based on the provided function.

    Parameters:
        points (np.ndarray): An array of points with shape (D, N).
        ham_func (Callable): The Hamiltonian function (either JIT-compiled or not).

    Returns:
        np.ndarray: An array of Hamiltonians with shape (D, 2, 2).

    Raises:
        RuntimeError: If both Numba and pure Python computations fail.
    """
    try:
        print(f"Attempting Numba-optimized computation for function '{ham_func.__name__}'")
        return numba_compute_hamiltonian(points, ham_func)
    except Exception as e_numba:
        print(f"Falling back to pure Python computation for function '{ham_func.__name__}'")

    try:
        return python_compute_hamiltonian(points, ham_func)
    except Exception as e_python:
        print(
            f"Error: Pure Python computation failed for function '{ham_func.__name__}' "
            f"with error: {e_python}"
        )
        raise RuntimeError(
            f"Hamiltonian computation failed for function '{ham_func.__name__}' using both "
            "Numba and Python implementations."
        )
