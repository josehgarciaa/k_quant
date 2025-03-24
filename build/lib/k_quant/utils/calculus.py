import numpy as np

def cumulative_integral(E, f_E):
    """
    Computes the cumulative integral of f(E) using the trapezoidal rule.

    Parameters:
        E (numpy.ndarray): Array of energy values (E_0, E_1, ..., E_n).
        f_E (numpy.ndarray): Array of function values f(E_0), f(E_1), ..., f(E_n).

    Returns:
        numpy.ndarray: Cumulative integral from E_0 to each E_i.
    """
    # Compute differences between consecutive energy points
    dE = np.diff(E)
    
    # Compute trapezoidal areas for each interval
    trap_areas = (f_E[:-1] + f_E[1:]) / 2 * dE
    
    # Cumulative sum to get the integral up to each E_i
    cumulative_integral = np.concatenate(([0], np.cumsum(trap_areas)))

    return cumulative_integral