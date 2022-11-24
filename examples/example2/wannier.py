import numpy as np
import multiprocessing as mp


class WeightedEdge:
    
    kpoint   = (0.0,0.0,0.0);
    direction= (0,0,0)
    weight   = 0.0j;
   
    def __init__(self, direction,weight):
        self.direction = direction;
        self.weight    = weight;
    
    def kpoint(self, kpoint):
        self.kpoint = kpoint;

    def bloch_hopping(self):
        self.kpoint = self.weight*np.exp(-2j*np.pi*np.dot( self.kpoint, self.direction) );

   
        




class WannierSystem:
    """A class that handles Wannier-based Tight binding models in different Wannier 


    :Keyword Arguments:
        The arguments are presented in order of preference, i.e, wout will be used over label. 
        
        label (str): The label of a Wannier90 calculation. This label will be used to look for the label.uc, label_hr.dat and label.xyz files. The positions and primitive lattice vectors in
        label.xyz and label.uc respectively, should be in the same units
        
        Convention
        It is important to highlight that in our label.xyz we only include the positions of the Wannier center and not the atoms
        
        
        """

    label = None;
    orbital_positions = None;
    primitive_vectors = None;
    wann_hamiltonian =None;

     
    def __init__(self, **kwargs ):

        #When the user submit a Wannier90 input, use it to initialice the system
        if "label" in kwargs:
            self.label= kwargs["label"];
            self.primitive_vectors = self.PrimitiveLatticeVectors(self.label+".uc");
            self.orbital_positions = self.OrbitalPositions(self.label+".xyz");
            self.wann_hamiltonian  = self.WannierHamiltonian(self.label+"_hr.dat");
            
            #Transform the positions in fractional units
            to_frac = np.linalg.inv(self.primitive_vectors).T;
            for k,v in self.orbital_positions.items():
                self.orbital_positions[k] = np.dot(to_frac,v);
            

    def OrbitalNumber(self):
        """Returns the number of orbitals in the model

        :return: Total number of orbitals in the model
        :rtype: int
        """
        return len(self.orbital_positions);


    def WannierHamiltonian(self,filename):
        """Load the wannier hamiltonian from a _hr.dat file

        :param filename: A string containing the path to the _hr.dat file to be loaded
        :type filename: str
        :return: a tuple containing the unit cell shifts, the row and c
        :rtype: _type_
        """

        file = open(filename, "r");
        description = file.readline();
        num_orbs = int(file.readline());
        numKPT  = int(file.readline());
        numKPT_lines = int(np.ceil(numKPT/15));
        KPT =[];
        for n in range(numKPT_lines):
            KPT.append(file.readline());
        file.close();

        assert num_orbs == self.OrbitalNumber(), "In WannierHamilton function the number of orbitals from the file do not coincide with the number of orbitals in the model " ;

        #Read and format automatically the data
        wann_data= np.genfromtxt(filename, skip_header= numKPT_lines+3)   


        #For proper packaging, it is efficient to know in advance the maximum number of hoppings
        #An efficient way to do this is to use a set. We will call this object, edges, from
        #its paralellism with the edges in a graph. 
        edges    = set( map(tuple,wann_data[:,0:3] ) );
        edge2idx = { e:i for i,e in enumerate(edges) };
        num_edges= len(edges);
        wann_ham = np.full( (num_orbs,num_orbs,num_edges), WeightedEdge((0.,0.,0.), 0.j) )  ;

        for s1,s2,s3,i,j,rv,iv in wann_data:
            edge = (s1,s2,s3);
            i,j  = int(i)-1, int(j)-1 ;
            
            assert edge in edge2idx, "Trying to include a hopping with an edge not present in the set of edges"           
            wedge = WeightedEdge(edge, rv+1j*iv);
            wann_ham[i,j, edge2idx[edge] ] = wedge;           
             
        return wann_ham;


    def OrbitalPositions(self, filename):
        """Load the orbital position from a xyz file 

        :param filename : A string containing the path to the xyz file to be loaded  
        :type filename: str
        """        
        file = open(filename, "r");
        numOrbs = int(file.readline());
        description = file.readline();

        xyz_data = list();
        for line in file:
            line = line.replace(" ","\t").replace("\n","").split("\t") 
            xyz_data.append([ elem for elem in line if elem!=""])
        
        orbital_positions = { oid+"_id"+str(i):(float(x),float(y),float(z)) for i, (oid, x,y,z) in enumerate(xyz_data) }
        assert numOrbs == len(orbital_positions), "In OrbitalPosition function the number of orbitals from the file is not compatible with the number of positions " ;

        return orbital_positions;         


    def PrimitiveLatticeVectors(self, filename):
        """Load the primitive lattice vectors from a win file
        :param filename : A string containing the path to the uc file to be loaded  
        :type filename: str

        """
        
        return np.loadtxt(filename)


    def hamiltonian(self,k):
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

        print( self.wann_hamiltonian[0,0,:] )



        #numOrbs = self.OrbitalNumber()
        #kpoints = np.array(kpoints, ndmin=1)
        #if kpoints.ndim == 1:
        #    kpoints = kpoints.reshape((1, -1))

        #H = np.zeros((kpoints.shape[0], numOrbs, numOrbs), dtype=complex)
        #tmp_array = np.empty_like(H)
        #for R, hop in self.wann_hamiltonian.items():
            # When the hopping matrices are very large, allocating new
            # arrays for the result of this multiplication (which is
            # of size len(k_array) * self.size**2) becomes expensive.
            # To avoid this, we reuse the same temporary array - even
            # if this is _slightly_ slower for single k-point calculations.
        #    np.multiply(
        #        np.exp(2j * np.pi * np.dot(kpoints, R)).reshape((-1, 1, 1)),
        #        self._array_cast(hop)[np.newaxis, :, :],
        #        out=tmp_array,
        #    )
        #    H += tmp_array
        #H += H.conjugate().transpose((0, 2, 1))
        #pos_exponential = np.array(
        #        [[np.exp(2j * np.pi * np.dot(kpoints, p)) for p in self.orbital_position]] ).transpose((2, 0, 1))
        #H = pos_exponential.conjugate().transpose((0, 2, 1)) * H * pos_exponential

