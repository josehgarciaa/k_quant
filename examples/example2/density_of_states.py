import numpy as np

import wannier as wan


wannsyst = wan.WannierSystem(label="linear_chain")


wannsyst.Hamiltonian(k = [1,20,300])