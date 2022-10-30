import numpy as np
from numpy.linalg import eigh, eigvalsh, norm

class bandstructure:
    
       
    def __init__(self, lat_vec, H_k ):
              
        self.lat_vec   = lat_vec;
        self.bandpath  = None;
        self.H_k       = H_k;

    def set_hamiltonian_k(self, H_k):
        self.H_k = H_k;

    def hamiltonian_k(self, k):
        return self.H_k(k);

    def set_bandpath(self,bandpath):
        self.bandpath = bandpath;

    def bandpath(self):
        return self.bandpath;


    def Momentum_Rec2AbsMatrix(self ):
        return 2*np.pi* np.linalg.inv(self.lat_vec);
    
    def toAbsoluteCoords(self,x):
             return np.dot( x, np.transpose(self.Momentum_Rec2AbsMatrix() ) );
                

    def band_kpoints(self , absolute_coords = True): 
        labels, segments, npoints = list(zip(*(self.bandpath)));
        segments = np.array(segments);
        last_pt  = segments[-1]; 
        segments = np.transpose( [ segments[:-1],segments[1:] ], axes=(1,0,2) )

        kpoints  = [ np.linspace( *seg, n,endpoint=False ) for seg,n in zip(segments,npoints) ] ;
        kpoints  = np.concatenate( [*kpoints,[last_pt]] , axis=0)
        #If required, rescale to absolute value
        if absolute_coords  is True :
            kpoints = np.dot( kpoints, np.transpose(self.Momentum_Rec2AbsMatrix() ) );

        return np.array(kpoints);

    def XLabels(self ):
        Xaxis = self.Xaxis()    
        labels, segments, npoints = list(zip(*(self.bandpath)));
        xpos = np.concatenate([[0], np.cumsum( npoints)])[:-1]

        return [list(Xaxis[xpos]),labels]                
       
    def Xaxis(self ):
        Xaxis=np.cumsum(np.linalg.norm(np.diff( self.band_kpoints() ,axis=0, prepend = 0),axis=1));
        Xaxis-=Xaxis[0];#Remove the initial value.
       
        return Xaxis;
    

    def operator_k(self, operator, kpoint ):
        w, v = np.linalg.eigh( self.H_k(kpoint) );

        operator = np.array(operator(kpoint) if callable(operator) else operator, dtype=complex);
        try:
            ops_k = np.sum(np.conj(v)* operator.dot(v), axis=0) ;
            return np.real(ops_k);
        except:
            print("problem with operators_k");
                
        return None;

    
    def compute_bands(self, fermi_energy = 0.0, proj_ops = None ):
        if proj_ops is None:
            bands = np.transpose( [ self.operator_k( self.H_k(kp), kp)   for kp in  self.band_kpoints() ]) - fermi_energy;
            return bands;

        ops = [self.H_k, *proj_ops];
        proj_bands = np.array( [ [ self.operator_k(op, kp) for op in  ops] for kp in  self.band_kpoints()] );
        
        return np.transpose(proj_bands, axes=(2,1,0));

    def energies_inwindow(bands, energy_window):
        Emin,Emax = energy_window;
        return np.any(bands>Emin,axis=1)*np.any(bands<Emax,axis=1);  

