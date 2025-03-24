from k_quant.utils.space import create_mesh
import numpy as np
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



def set_lattice_vector( __lat_vec):
    set_param("global_lattice_vectors", __lat_vec)

def get_lattice_vector( ):
    return get_param("global_lattice_vectors")


def set_kmesh( dims):
    
    mesh = create_mesh(dims)
    global_lat = get_lattice_vector()
    rec2cart = 2*np.pi* np.linalg.inv(global_lat).T
    kmesh =np.dot( mesh, rec2cart )
    print(kmesh.dot(global_lat[1]))
    set_param("kmesh", np.dot( mesh, rec2cart ) )

    
def get_kmesh():
    return get_param("kmesh", default=None)

    
def set_energy_grid(energies):
    set_param("energy_grid", energies)
    print("Create a grid")
    
def get_energy_grid():
    return get_param("energy_grid", default=None)
