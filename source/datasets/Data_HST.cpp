#include <imcmc/parser++.hpp>
#include "HST.hpp"

using namespace std;
using namespace imcmc::parser;

Data_HST::Data_HST() {
//  default values
    H0          = 73.8;
    sigma_H0    = 2.4;
}

Data_HST::~Data_HST() {
// do nothing
}


void Data_HST::Init(std::string& dataset) {

    data_info.GetInfo(dataset);
    H0 = Read::Read_Double_from_File(dataset, "H0");
    sigma_H0= Read::Read_Double_from_File(dataset, "sigma_H0");

}
