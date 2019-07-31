#ifndef __MISC__
#define __MISC__

#include <iostream>
#include <string>
#include <map>

void MPI_cout( std::string mesg, int rank=0 );

// // typedef std::map<std::string, double>  imcmc_double;

void Print_ImcmcDouble(std::map<std::string, double>& param);


#define ClassMC_ERROR(ERR_MSG) { std::cout 		\
<< "\t#### ClassMC Error ####\n"				\
<< "#--- File Name: " << __FILE__ << "\n" 		\
<< "#--- Line    #: " << __LINE__ << "\n" 		\
<< "#--- Func Name: " << __FUNCTION__ << "\n" 	\
<< "#--- Error Msg: " << ERR_MSG << "\n\n"; 	\
exit(0);										\
}

#define ClassMC_WARNING(WARNING_MSG) { std::cout\
<< "\t#### ClassMC Warning ####\n"				\
<< "#--- File Name: " << __FILE__ << "\n" 		\
<< "#--- Line    #: " << __LINE__ << "\n" 		\
<< "#--- Func Name: " << __FUNCTION__ << "\n" 	\
<< "#--- Warning Msg: " << WARNING_MSG << "\n\n"; 	\
}


#define MPI_ClassMC_ERROR(ERR_MSG) {        			\
	if( MPI::COMM_WORLD.Get_rank() == 0 ){              \
		std::cout 										\
		<< "\t#### ClassMC Error ####\n"				\
		<< "#--- File Name: " << __FILE__ << "\n" 		\
		<< "#--- Line    #: " << __LINE__ << "\n" 		\
		<< "#--- Func Name: " << __FUNCTION__ << "\n" 	\
		<< "#--- Error Msg: " << ERR_MSG << "\n\n"; 	\
	}													\
	MPI::Finalize();								\
	exit(0);										\
}

#define MPI_ClassMC_WARNING(WARNING_MSG) {     			\
	if( MPI::COMM_WORLD.Get_rank() == 0 ){              \
		std::cout 										\
		<< "\t#### ClassMC Warning ####\n"				\
		<< "#--- File Name: " << __FILE__ << "\n" 		\
		<< "#--- Line    #: " << __LINE__ << "\n" 		\
		<< "#--- Func Name: " << __FUNCTION__ << "\n" 	\
		<< "#--- Warning Msg: " << WARNING_MSG << "\n\n"; 	\
	}														\
}


#define DEBUG_STOP {						\
	std::cout << "\t------------------";	\
	std::cout << "\t> MPI debug stop <\n"; 	\
	std::cout << "\t------------------";	\
	exit(0);								\
}

#define MPI_DEBUG_STOP {						\
	if( MPI::COMM_WORLD.Get_rank() == 0 ){      \
		std::cout << "----------------------------------------\n";	\
		std::cout << "> =========== MPI debug stop ==========<\n"; 	\
		std::cout << "> File Name: " << __FILE__ << "\n";             \
		std::cout << "> Line    #: " << __LINE__ << "\n";             \
		std::cout << "> Func Name: " << __FUNCTION__ << "\n";         \
		std::cout << "----------------------------------------\n";	\
	}											\
	MPI::Finalize();							\
	exit(0);									\
}


//	Macros defined to count time (in seconds)
#define DEBUG_TIME_START(_time_) {  \
	std::cout << "TIMER: start counting ...\n";\
    _time_ = std::time(NULL);  		\
}

#define DEBUG_TIME_END(_time_) { 				\
    std::cout   								\
    << "=====================================\n"\
    << "DEBUG: Wall time passed: " 				\
    << std::difftime(std::time(NULL), _time_) 	\
    << " seconds.\n" 							\
    << "=====================================\n";\
}

#endif // __MISC__

