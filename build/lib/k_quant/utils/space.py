import numpy as np    
    
    
    
    
    
def create_mesh(dims) :
    meshkgrid = np.meshgrid(*[np.linspace(0, 1, d, endpoint=False) for d in dims], indexing='ij')
    return np.transpose([x.flatten() for x in meshkgrid])