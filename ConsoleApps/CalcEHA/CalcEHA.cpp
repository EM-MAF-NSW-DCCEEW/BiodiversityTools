/*
CalcEHA.cpp - Console app for running a context CBA analysis to produce EHA
Copyright(C) 2024 State of New South Wales and Department of Climate Change, Energy, the Environment and Water (DCCEEW)
Author: Jamie Love, Ecological Modelling, Science and Insights Division

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <https://www.gnu.org/licenses/>.
*/

#include <fstream>
#include <iostream>

#include "FileUtilsLib.h"
#include "SpatialContextLib.h"

//TODO rename project EHA_CBA?
int main(int argc, char *argv[])
{
	//Options
	//Pause after use
	if (argc < 2) {
		std::cout << "\nUsage: " << argv[0] << " ParameterFile.par \n"
			<< "\nIf ParameterFile.par exists it is processed, otherwise it is created as a template.\n\n";
		return -1;
	}
	std::fstream testFS;
	testFS.open(argv[1], std::ios::in);
	    if (!testFS.is_open()) {
		std::cout << "\nThe parameter file " << argv[1] << " does not exist and is being created as a template\n\n";
		if (!SetParam(argv[1], "INPUTS", "HAB_FN", "Filename. Habitat input grid file name without extension.") || 
			!SetParam(argv[1], "INPUTS", "PRM_FN", "Filename. Permeability input grid file name without extension.") || 
		    !SetParam(argv[1], "OUTPUTS", "EHA_FN", "Filename. EHA output filename without extension.") || 
		    !SetParam(argv[1], "WINDOW", "LOOKUP_FN", "Filename. Petals look up file name with extension. Source and destination window sizes are derived from this.") || 
		    !SetParam(argv[1], "WINDOW", "TRANSLATE_FN", "Filename. Petals translate file name with extension. Source and destination window sizes are derived from this.") || 
		    !SetParam(argv[1], "ANALYSIS", "MULTFOCAL", "true/false, Whether to multiply connectivity by the focal cell habitat value. Optional, default = false.") || 
		    !SetParam(argv[1], "ANALYSIS", "FOCALPOWER", "Floating point number. Exponent applied to the focal cell habitat value. Optional, default = 1.0.") || 
		    !SetParam(argv[1], "ANALYSIS", "SUMPOWER", "Floating point number. Exponent applied to occupancy and MPC. Optional, default = 1.0.")) 
			return msgErrorStr(-99, "Unable to create parameter file template:", argv[1]);
		return -2;
	}
	testFS.close();
	int returnValue = EHA_CBA(argv[1]);
	if (returnValue != 0) std::cout << "Error code: " << returnValue << "\n";
	else std::cout << "Complete!" << "\n";
	return returnValue;
}