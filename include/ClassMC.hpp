
//	=====================
//	standard C++ headers
//	=====================

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

//	=====================
//	OpenMPI header
//	=====================

#ifndef __CLASSMC_MPI__
#define __CLASSMC_MPI__
	#include "mpi.h"
#endif

//	=====================
//	ClassEngine headers
//	=====================

#include "class.h"
#include "Engine.hh"
#include "ClassEngine.hh"


//	========================================
//	Model, Data, ParamNames and Likelihoods
//	========================================

#include "CosmoTheory.hpp"
#include "DataList.hpp"
#include "ParamList.hpp"
#include "Likelihoods.hpp"

//	================================
//	Armadillo C++ linear algebra lib
//	================================
#include <armadillo>