/* 
 Copyright 2018 Altobelli de Brito Mantuan.
 Distributed under the GNU General Public License, Version 3.
 Last update : August 25th, 2018.
*/

#include <iostream>
#include <chrono>

#undef NDEBUG
#include <assert.h>

#define EIGEN_RUNTIME_NO_MALLOC // Define this symbol to enable runtime tests for allocations

#include "ECMI.h"

int main(int argc, char* argv[]){

    if (argc < 4){
        std::cout << "Parameter missing. Please follow the instructions below: " << std::endl;
		std::cout << "SCIM.exe [INPUT_PATH] [DR] [OUTPUT_PATH]" << std::endl;
        return 1;
    }

    double prc_corr = 1 - atof(argv[2]); // dr = 0.0 the automatic radio region and dr 1.0 = bigest radio region of cluster
	std::string strFile = argv[1];
    std::string strFileOut = argv[3];

    std::cout << "file = " << strFile << std::endl;
    std::cout << "parameter correlation precision= " << atof(argv[3]) << std::endl;
	std::cout << "save result in = " << strFileOut << std::endl;

    ECMI ecmi(strFile, strFileOut);
    ecmi.runAlgorithm(prc_corr);
	std::cout << std::endl;
	std::cout << "Total itemset = " << ecmi.getQtdItemsetGerados() << std::endl;
    
    return 0;
}
