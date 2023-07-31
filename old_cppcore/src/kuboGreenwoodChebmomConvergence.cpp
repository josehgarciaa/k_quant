// C & C++ libraries
#include <iostream> /* for std::cout mostly */
#include <string>   /* for std::string class */
#include <fstream>  /* for std::ofstream and std::ifstream functions classes*/
#include <stdlib.h>
#include <chrono>


#include "kpm_noneqop.hpp" //Message functions
#include "chebyshev_moments.hpp"
#include "sparse_matrix.hpp"
#include "quantum_states.hpp"
#include "chebyshev_solver.hpp"

int main(int argc, char *argv[])
{
	if ( !(argc == 7) )
	{
		chebyshev::convergence::printHelpMessage();
		return 0;
	}
	else
		chebyshev::convergence::printWelcomeMessage();
	
	const std::string
		LABEL = argv[1],
		S_OPR = argv[2],
		S_OPL = argv[3],
		S_NMOM= argv[4],
		S_ETA = argv[5],
		S_E0  = argv[6];

	const int numMoms  = atoi(S_NMOM.c_str() );
	const double eta   = stod( S_ETA );
	const double energ = stod( S_E0 );
	chebyshev::Moments1D chebMoms(numMoms); //load number of moments


	SparseMatrixType OP[3];
	OP[0].SetID("HAM");
	OP[1].SetID(S_OPR);
	OP[2].SetID(S_OPL);

	// Build the operators from Files
	SparseMatrixBuilder builder;
	std::array<double,2> spectral_bounds;	
	for (int i = 0; i < 3; i++)
	{
		std::string input = "operators/" + LABEL + "." + OP[i].ID() + ".CSR";
		builder.setSparseMatrix(&OP[i]);
		builder.BuildOPFromCSRFile(input);
	
		if( i == 0 ) //is hamiltonian
		//Obtain automatically the energy bounds
		 spectral_bounds = chebyshev::utility::SpectralBounds(OP[0]);
	};
	//CONFIGURE THE CHEBYSHEV MOMENTS
	chebMoms.SystemLabel(LABEL);
	chebMoms.BandWidth ( (spectral_bounds[1]-spectral_bounds[0])*1.0);
	chebMoms.BandCenter( (spectral_bounds[1]+spectral_bounds[0])*0.5);
	chebMoms.SystemSize(OP[0].rank() );
	chebMoms.SetHamiltonian(OP[0]);

	chebMoms.Print();

	//Define thes states youll use
	//Factory state_factory ;
	
	std::string outputfilename="ConvNonEqOp"+S_OPR+"-"+S_OPL+LABEL+"KPM_M"+S_NMOM+"eta"+S_ETA+"E"+S_E0+".dat";
	chebyshev::sequential::KuboGreenwoodChebMomConvergence( energ,eta,OP[1], OP[2], chebMoms );

	std::cout<<"Saving convergence data in "<<outputfilename<<std::endl;
	std::ofstream outputfile( outputfilename.c_str() );
	for( int m =0; m < chebMoms.HighestMomentNumber(); m++ )
		outputfile << chebMoms(m).real() <<std::endl;
	outputfile.close();


	std::cout<<"End of program"<<std::endl;
	return 0;
}

