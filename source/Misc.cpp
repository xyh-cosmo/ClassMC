#include "Misc.hpp"

void MPI_cout( std::string mesg, int rank ){
    if( rank == 0 ){
        std::cout << "==> " << mesg << std::endl;
    }
}


// typedef std::map<std::string, double>  imcmc_double;

void Print_ImcmcDouble(std::map<std::string, double>& param){

	std::map<std::string, double>::iterator it;

	for( it = param.begin(); it != param.end(); ++it ){
		std::cout << "debug ==> " << it->first << " = " << it->second << std::endl;
	}
}
