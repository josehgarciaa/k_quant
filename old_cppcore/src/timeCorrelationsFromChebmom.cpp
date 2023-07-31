
// C & C++ libraries
#include <iostream>		/* for std::cout mostly */
#include <string>		/* for std::string class */
#include <fstream>		/* for std::ofstream and std::ifstream functions classes*/
#include <vector>		/* for std::vector mostly class*/
#include <complex>		/* for std::vector mostly class*/
#include <omp.h>

#include "chebyshev_solver.hpp"
#include "chebyshev_coefficients.hpp"
#include "chebyshev_moments.hpp"

void printHelpMessage();
void printWelcomeMessage();

int main(int argc, char *argv[])
{	
	if (argc != 4)
	{
		printHelpMessage();
		return 0;
	}
	else
		printWelcomeMessage();

	const double broadening  = stod(argv[2]);
	const double fermiEner  = stod(argv[3]);


	//Read the chebyshev moments from the file
 	chebyshev::MomentsTD mu(argv[1]); 
	//and apply the appropiated kernel
	mu.ApplyJacksonKernel(broadening);
	
	mu.Print();
	

	const int num_div = 30*mu.HighestMomentNumber();

	const double
	x = mu.Rescale2ChebyshevDomain(fermiEner);
		
	
	std::cout<<"Computing the equillibriu, time-correlations using "<<mu.HighestMomentNumber()<<" x "<<mu.MaxTimeStep()<<" TD moments "<<std::endl;
	std::cout<<"on the fermi energy:"<<fermiEner<<" which is normalized to "<<x<< std::endl;
	std::cout<<"The first 10 moments are"<<std::endl;
	for(int m0 =0 ; m0< 1; m0++ )
	for(int m1 =0 ; m1< 10; m1++ )
		std::cout<<mu(m0,m1).real()<<" "<<mu(m0,m1).imag()<<std::endl;
	
	std::string
	outputName  ="mean"+mu.SystemLabel()+"EF"+std::to_string(fermiEner)+"JACKSON.dat";


	std::cout<<"Saving the data in "<<outputName<<std::endl;
	std::cout<<"PARAMETERS: "<< mu.SystemSize()<<" "<<mu.HalfWidth()<<std::endl;

	std::ofstream outputfile( outputName.c_str() );
	for( int n = 0 ; n < mu.MaxTimeStep(); n++)
	{
		
//		if( n == mu.MaxTimeStep()-1)for( double x=-0.99; x<0.99; x=x+0.01){
		double output = 0.0;
		for( int m = 0 ; m < mu.HighestMomentNumber() ; m++)
				output += delta_chebF(x,m)*mu(m,n).real() ;
		output *=  mu.SystemSize()/mu.HalfWidth();
//		outputfile<<x*mu.HalfWidth() + mu.BandCenter() <<" "<<output <<std::endl;}
		outputfile<<n*mu.TimeDiff() <<" "<<output <<std::endl;
		
	}
	outputfile.close();

	std::cout<<"The program finished succesfully."<<std::endl;
return 0;
}
	

void printHelpMessage()
{
	std::cout << "The program should be called with the following options: moments_filename broadening(meV) maximum_time (fs)" << std::endl
			  << std::endl;
	std::cout << "moments_filename will be used to look for .chebmomTD file" << std::endl;
	std::cout << "broadening in (meV) will define the broadening of the delta functions" << std::endl;
	std::cout << "The Fermi energy in (meV) ." << std::endl;

};

void printWelcomeMessage()
{
	std::cout << "WELCOME: This program will sum the chebyshev moments for obtaining the time-evolution of a correlator in equilibrium" << std::endl;
};
