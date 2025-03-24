


def validate_operator_dimensions(operators):
    """
    Check that all matrices in the list have the same dimensions.

    Parameters:
    -----------
    operators: list of np.ndarray
        A list of matrices to check.

    Returns:
    --------
    tuple:
        The shape of the matrices if all are identical.

    Raises:
    -------
    ValueError:
        If matrices do not all have the same shape.
    """
    if not operators:
        raise ValueError("The list of operators is empty.")

    shape = operators[0].shape
    for op in operators[1:]:
        if op.shape != shape:
            raise ValueError(f"Operator shapes do not match. Expected {shape}, but got {op.shape}.")

    return shape