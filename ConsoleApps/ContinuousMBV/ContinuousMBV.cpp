// ContinuousMBV.cpp : Defines the entry point for the console application.
//

#include <fstream>
#include <iostream>

#include "FileUtilsLib.h"
#include "SpatialContextLib.h"

int main(int argc, char *argv[])
{
	//Options
	//Pause after use
	if (argc < 2) {
		std::cout << "\nUsage: " << argv[0] << " ParameterFile.par \n"
			<< "\nIf ParameterFile.par exists it is processed, otherwise it is created as a template.\n\n";
		return -1;
	}
	std::string argv1(argv[1]);
	std::fstream testFS;
	testFS.open(argv[1], std::ios::in);
	if (!testFS.is_open()) {
		std::cout << "\nThe parameter file " << argv[1] << " does not exist and is being created as a template\n\n";
		if (!SetParam(argv1, "OUTPUTS", "OUT_FN", "Filename. Resilience denominator output filename without extension. ") ||
			!SetParam(argv1, "RESILIENCE", "HAB_FN", "Habitat grid used in the numerator") ||
			!SetParam(argv1, "RESILIENCE", "NTXGRIDS", "Integer. The number of Transform grid pairs TXnSRC, TXnDST where n in 0..NTXGRIDS-1.") ||
			!SetParam(argv1, "RESILIENCE", "MAXSAMPLES", "Integer. The maximum number of pixels to sample from TXDST. Optional, default use SAMPLERATE.") ||
			!SetParam(argv1, "RESILIENCE", "SAMPLERATE", "Integer. The rate at which pixels are sampled from TXDST (1 in n). Optional, default use MAXSAMPLES.") ||
			!SetParam(argv1, "RESILIENCE", "TX0SRC", "Filename. The first transform grid used at each cell (i). Additional TXnSRC required for n in 0..NTXGRIDS-1.") ||
			!SetParam(argv1, "RESILIENCE", "TX0DST", "Filename. The first transform grid used at each sample (j). Additional TXnDST required for n in 0..NTXGRIDS-1."))
			return msgErrorStr(-99, "Unable to create parameter file template:", argv[1]);
		return -2;
	}
	testFS.close();
	int returnValue = ContinuousMBV(argv[1]);
	if (returnValue != 0) std::cout << "Error code: " << returnValue << "\n";
	else std::cout << "Complete!" << "\n";
	return returnValue;
}