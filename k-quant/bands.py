import numpy as np
import numpy.typing as npt
from numpy.linalg import eigh, eigvalsh, norm
from numpy import exp, cos, sin , pi, kron


class bands:
    """This class will compute the bandstructure for different inputs of models

    The __init__ method may be documented in either the class level
    docstring, or as a docstring on the __init__ method itself.

    Either form is acceptable, but the two should not be mixed. Choose one
    convention to document the __init__ method and be consistent with it.

    Note:
        Do not include the `self` parameter in the ``Args`` section.

    Args:
        msg (str): Human readable string describing the exception.
        code (:obj:`int`, optional): Error code.

    Attributes:
        msg (str): Human readable string describing the exception.
        code (int): Exception error code.

    """
       
    def __init__(self, label ):
              
        self.lat_vec   = None;
        self.bandpath  = np.array( ['G',1,[0,0,0]] );
        self.H_k       = None;

    def set_hamiltonian_k(self, H_k):
        """_summary_

        Args:
            H_k (_type_): _description_
        """
        self.H_k = H_k;

    def hamiltonian_k(self, k):
        """_summary_

        Args:
            k (_type_): _description_

        Returns:
            _type_: _description_
        """
        return self.H_k(k);

    def set_bandpath(self,bandpath):
        """_summary_

        Args:
            bandpath (_type_): _description_
            absolute_coords (bool, optional): _description_. Defaults to False.
        """
        self.bandpath = bandpath;

    def bandpath(self):
        """_summary_

        Args:
            bandpath (_type_): _description_
            absolute_coords (bool, optional): _description_. Defaults to False.
        """
        return self.bandpath;


    def Momentum_Rec2AbsMatrix(self ):
        return 2*np.pi* np.linalg.inv(self.lat_vec);
    
    #If required, rescale to absolute value
    def toAbsoluteCoords(self,x):
             return np.dot( x, np.transpose(self.Momentum_Rec2AbsMatrix() ) );
                

    def band_kpoints(self , absolute_coords = False): 
        """Compute the k-points defined in the bandpath.  

        Args:
            absolute_coords (bool, optional): Determine whereas the k-points are given in cartesian (in units of 2pi/a) or reciprocal units. Defaults to False.

        Returns:
            ndarray: An array of k-points
        """
        kpoints   = list(); 
        init_k = self.bandpath[0][2];
        for path_label, npoint, end_k in self.bandpath[1:]: #Not consider initial point anymore
                path_kpoints = np.linspace(init_k,end_k,npoint, endpoint=False ); 
                init_k  = end_k
                for kp in path_kpoints:
                    kpoints.append(kp)
        kpoints.append(init_k);

        #If required, rescale to absolute value
        if absolute_coords  is True :
             kpoints = np.dot( kpoints, np.transpose(self.Momentum_Rec2AbsMatrix() ) );
                
        return np.array(kpoints);


    def XLabels(self ):
        """Get the x-axis labels associated with the bandpath

        Returns:
            list: _description_
        """
        Xaxis = self.bandsXaxis()
        xpos = (np.cumsum(self.bandpath[:,1],dtype=int)-1)
        return (Xaxis[xpos],self.bandpath[:,0])                
    
    
    def bandsXaxis(self ): #The kpoints are assume to be in absolute coordinates 
        kpoints=self.band_kpoints(absolute_coords = True);
        Xaxis=np.cumsum(np.linalg.norm(np.diff( kpoints,axis=0, prepend = 0),axis=1));
        Xaxis-=Xaxis[0];#Remove the initial value.
        return Xaxis;
    
    
    def compute_band_structure(self, fermi_energy = 0.0, proj_op = None, ax=None, plot_proj=False, proj_range=None ):
     
        #Compute the kpoints based on the path
        kpoints= self.band_kpoints();

        #Compute the eigenvalues and the projected values 
        peigenvals  = self.compute_dispersion( kpoints , proj_op)  ;

        if proj_op is None:
            bands = peigenvals.T- fermi_energy;
            return self.plot_band_structure( bands = bands, ax=ax );

        bands = peigenvals[:,0].T - fermi_energy;
        projs = peigenvals[:,1].T;
        
        return self.plot_band_structure( bands = bands , projs = projs , ax=ax, plot_proj = plot_proj, proj_range=proj_range);

    
    def plot_band_structure(self, bands, projs = None, ax=None, plot_proj = False, proj_range=None):
        
        xaxis  = self.bandsXaxis();
            
        #PLOTING
        #plot options
        if  ax is None:
            fig = plt.figure();
            _ax = fig.add_subplot();
        else:
            _ax  = ax;

        _ax.tick_params(axis='both', which='major', labelsize=16);
        
        if plot_proj is True:
            for i,proj in enumerate(projs):
                _ax.plot(xaxis,proj );
        else:
            _ax.set_ylabel("Energy (eV)", fontsize=16);
            vmin,vmax = 0,1;
            if projs is not None:
                vmin = np.min(projs);
                vmax = np.max(projs);
                if proj_range is not None:
                    vmin,vmax = list(proj_range);
            for i,band in enumerate(bands):
                _ax.plot(xaxis,band , color = 'k', lw=band_lw);
                c = 'k';
                if projs is not None:
                    c = projs[i];
                #Add points to the plot
                im = _ax.scatter(xaxis, band,s=proj_s, c=c,cmap=cmap_name,vmin=vmin, vmax=vmax);
                im.set_facecolor("none");        
            #Create falso plot for color bar
            if projs is not None and ax is None: #Add the color bar whenever you have the a projectedprlt
                im = _ax.scatter(xaxis, bands[0],s=0, c=np.linspace( vmin,vmax,len(bands[0]) ),cmap=cmap_name);
                fig.colorbar(im, ax=_ax); 

        klabels, path_labels = self.get_XLabels()
        _ax.set_xticks(klabels);
        _ax.set_xticklabels(path_labels);
                
        if ax is None:

            return fig,_ax;

        ax = _ax;
        return ax;

    def compute_dispersion(self,kpoints, proj_op=None, absolute_coords = False):
        if ( absolute_coords is True ):
            Abs2Rec = np.linalg.inv(np.transpose(self.Momentum_Rec2AbsMatrix() ) );
            kpoints = np.dot( kpoints , Abs2Rec);
        
        return np.array( [self.projected_eigenvalues( kp, proj_op ) for kp in kpoints ] );        
      
    def compute_fermi_surface(self,kpoints, fermi_energy = 0.0, tol = None ,  proj_op = None ):

        #Compute the eigenvalues and the projected values 
        peigenvals  = self.compute_dispersion(kpoints, proj_op=proj_op)
        eigenvalues = np.array(list(map(list,peigenvals[:,0])))
        eigenvalues-= fermi_energy;#Compute the eigenvalues

        #Determine the important kpoints
        relevant_kpoints = np.any(np.abs(eigenvalues - fermi_energy) < tol, axis=1);
        #Select the kpoints
        kpoints = kpoints[relevant_kpoints];

        if proj_op is None:
            return kpoints,eigvals[relevant_kpoints];

        #If one requires a projection, select the kpoints and then use them to compute the projections
        
        peigvals = np.array(list(map(list,peigenvals[:,1])))
        return kpoints,peigvals[relevant_kpoints];

    
    def refine_kpoints(self, kpoints):

        nkp= 2*len(kpoints);
        kmin = np.array([ np.min(kpoints[:,0]),np.min(kpoints[:,1]),0 ] );
        kmax = np.array([ np.max(kpoints[:,0]),np.max(kpoints[:,1]),0 ] );
        kpoints= np.array(list(np.ndindex((nkp,nkp,1))))/(nkp-1)*(kmax-kmin) + kmin ;#kpoints in recpricola lattice vectors

        return kpoints;
    
    
    def compute_DOS(self, energies,gridpoints, fermi_energy = 0.0, broadening = 0.1,  kpoints = None ,proj_op = None ):
        
        if kpoints is None:
            n1,n2,n3 = gridpoints
            kpoints = ( np.mgrid[0:1:(n1+1)*1j, 0:1:(n2+1)*1j, 0:1:(n3+1)*1j].T)[:-1,:-1,:-1].reshape(n1*n2*n3,3);

        eta=broadening;
        def gaussian(x):
            return np.exp( -(x/eta)**2/2 )/eta/np.sqrt(2*np.pi)

        dim = len(kpoints);
        peigenvals  = self.compute_dispersion(kpoints, proj_op=proj_op)
        eigenvalues = np.array(list(map(list,peigenvals[:,0]))).flatten()
        eigenvalues -= fermi_energy;#Compute the eigenvalues
        
        
        if proj_op is None:
            return np.array( [np.sum(gaussian(eigenvalues-EF)) for EF in energies] )/dim;

        P = np.array(list(map(list,peigenvals[:,1]))).flatten()
        return np.array( [np.sum(P*gaussian(eigenvalues-EF)) for EF in energies] )/dim;


    def band_spin_texture(self, kpoints, bidx, uop = None, vop = None, zop=None, zlims = None, ax=None ):
        #Compute first the mean values using the model

        Ops= [uop, vop, zop ] ;

        def get_band_data( x, bidx ):
            y = None;
            if x.shape[0] == kpoints.shape[0]:
                y = x;
            else:
                y = self.compute_dispersion(kpoints, proj_op = x )[:,1];
            y = np.array(list(map(list,y))).T;
            return y[bidx];

        #If is a compatible array of data, add it directly to Ops
        Ops= [ get_band_data(op,bidx) for op in Ops ];

        uop, vop, zop = Ops;
        KX,KY,KZ = np.transpose(kpoints);

        # Create triangulation.
        triang = mtri.Triangulation(KX, KY)

        def get_mask( x, y, z ):
            return mtri.CubicTriInterpolator(triang, z, kind='geom')(x, y) > 0.;

        def xypoints( npts ):
            return np.meshgrid(np.linspace(KX.min(), KX.max(), npts), np.linspace(KY.min(), KY.max(), npts));

        # Set up the figure
        if( ax is None ):
            fig, ax = plt.subplots(dpi=400);

        # Interpolate to regularly-spaced quad grid.
        xi, yi = xypoints(npts=100)
        fzi    = mtri.CubicTriInterpolator(triang, zop, kind='geom');
        zi     = np.ma.array( fzi(xi,yi) , mask=get_mask(xi, yi,zop) )

        cs = ax.contourf(xi, yi, zi,levels=100, cmap="bone")
#        fig.colorbar(cs, ax=ax, shrink=0.9 , label=r"$\rm Energy (meV)$")
        # This is the fix for the white lines between contour levels
        for c in cs.collections:
            c.set_edgecolor("face")


        xi, yi = xypoints(npts=20)
        fui= mtri.CubicTriInterpolator(triang, uop, kind='geom');
        fvi= mtri.CubicTriInterpolator(triang, vop, kind='geom');
        ui = np.ma.array(fui(xi, yi), mask=get_mask(xi, yi,zop));
        vi = np.ma.array(fvi(xi, yi), mask=get_mask(xi, yi,zop));
        Q = ax.quiver(xi, yi, ui, vi, color="C1",  scale=15, width=0.022/4);

        return ax;
    
    
    def createRandomPhase(self, dims, output, mask=None):
        def linearize( i , dims ):
            x,y,z =i ;
            nx,ny,nz =dims ;
            return z*nx*ny + y*nx + x;

        nx,ny,nz = dims;
        numstates  = matdim = len(self.ham_operator([0,0,0]))
        dim = np.prod(dims);
        states = np.zeros((numstates,numstates*dim), dtype=complex)
        grid =(np.mgrid[0:nx,0:ny,0:nz].T).reshape( dim ,3);

        for kp in grid:
            if mask is not None and mask(kp/dims):
                continue;
            else:
                k = linearize(kp,dims);
                kp = kp/dims;
                w, vs  = np.linalg.eigh(self.ham_operator(kp))
                vs = vs.T;
                for i,v in enumerate(vs):
                    states[i,k*matdim:(k+1)*matdim]= np.exp(2j*np.pi*np.random.random())*v;

        for s in range(numstates):
            for i in range(matdim):
                S = states[s][i::matdim].reshape(dims).T
                FS = np.fft.fftn(S)
                states[s][i::matdim] = FS.T.flatten();
            states[s] /=np.linalg.norm(states[s]);

        f = open(output,'w')
        f.write("vector\n")
        f.write( str(numstates)+"\n")
        f.write( str(dim*matdim)+"\n")

        for state in states:
            for x in state:
                f.write(str(x.real)+" "+str(x.imag)+"\n")
        f.close()

    def projected_bands(self, proj_ops  ):
        #When no operator submited, return only the band structure
        if proj_op is None:
            return ( eigvalsh( Hk ) ) ;

        #If not, compute the eigen vectors
        w, v = np.linalg.eigh( Hk );#The column w[:, i] is the normalized eigenvector of v[i] eigenvalue

        w = np.diag( (np.conj(v.T) ) .dot(Hk.dot(v) ) ) ;
        p = np.diag( (np.conj(v.T) ) .dot(proj_op.dot(v) ) ) ;
        return ( np.real(w),np.real(p) );