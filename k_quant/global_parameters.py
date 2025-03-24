from k_quant.utils.space import create_mesh

_parameters = {}

def set_param(name, value):
    """Sets a global parameter."""
    _parameters[name] = value

def get_param(name, default=None):
    """Gets a global parameter. Returns default if not found."""
    return _parameters.get(name, default)

def reset_params():
    """Resets all global parameters."""
    _parameters.clear()

def list_params():
    """Lists all global parameters."""
    return dict(_parameters)


def set_kmesh( dims):
    set_param("kmesh", create_mesh(dims))

def get_kmesh():
    return get_param("kmesh", default=None)

    
def set_energy_grid(energies):
    set_param("energy_grid", energies)
    print("Create a grid")
    
def get_energy_grid():
    return get_param("energy_grid", default=None)
