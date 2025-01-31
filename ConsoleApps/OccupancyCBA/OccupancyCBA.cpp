/*
REMP_CBA.cpp - TODO Description
Copyright(C) 2024 State of New South Wales and Department of Climate Change, Energy, the Environment and Water (DCCEEW)
NSW DCCEEW Metrics and Forecasting Ecological Modelling Team 
ecological.modelling@environment.nsw.gov.au

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


// REMP_CBA.cpp : Defines the entry point for the console application.
//

//TODO standardise parameter naming convention

//#include "stdafx.h"
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
	std::fstream testFS;
	testFS.open(argv[1], std::ios::in);
	if (!testFS.is_open()) {
		std::cout << "\nThe parameter file " << argv[1] << " does not exist and is being created as a template\n\n";
		if (!SetParam(argv[1], "INPUTS", "HAB_FN0", "Filename. Habitat input grid file name without extension for time step 0. Additional _FN1..._FNt-1 required for t time steps.") || 
		    !SetParam(argv[1], "INPUTS", "PRM_FN0", "Filename. Permeability input grid file name without extension for time step 0. Additional _FN1..._FNt-1 required for t time steps.") || 
		    !SetParam(argv[1], "INPUTS", "NHA_FN0", "Filename. NHA input grid file name without extension for time step 0. Additional _FN1..._FNt-1 required for t time steps. Optional, default = use EXT_FN or E.") || 
		    !SetParam(argv[1], "INPUTS", "EXT_FN0", "Filename. Extinction input grid file name without extension for time step 0. Additional _FN1..._FNt-1 required for t time steps. Optional, default = use NHA_FN or E.") || 
		    !SetParam(argv[1], "INPUTS", "OCC_IN", "Filename. Initial occupancy factor input grid file name without extension. Optional, default = OCC_VAL or 1.0.") || 
		    !SetParam(argv[1], "INPUTS", "MPC_IN", "Filename. Initial MPC factor input grid file name without extension. Optional, default = MPC_VAL or 1.0.") || 
		    !SetParam(argv[1], "INPUTS", "OCC_VAL", "Floating point number. Initial occupancy factor value. Optional, default = OCC_IN or 1.0.") || 
		    !SetParam(argv[1], "INPUTS", "MPC_VAL", "Floating point number. Initial MPC factor value. Optional, default = MPC_IN or 1.0.") || 
		    !SetParam(argv[1], "WINDOW", "LOOKUP_FN", "Filename. Petals look up file name with extension. Source and destination window sizes are derived from this.") || 
		    !SetParam(argv[1], "WINDOW", "TRANSLATE_FN", "Filename. Petals translate file name with extension. Source and destination window sizes are derived from this.") || 
		    !SetParam(argv[1], "SPECIES", "LAMBDA", "Floating point number. lambda parameter. Optional, default derived from MCalcs.") || 
		    !SetParam(argv[1], "SPECIES", "C", "Floating point number. DEPRECIATED Use LAMBDA! Colonisation parameter c. Optional, default derived from MCalcs.") || 
		    !SetParam(argv[1], "SPECIES", "E", "Floating point number. DEPRECIATED Use LAMBDA! Extinction parameter e. Optional, default = 1.0.") || 
		    !SetParam(argv[1], "ANALYSIS", "TIMESTEPS", "Integer. Number of time steps to be processed. Must provide _FN1...FNt-1 Filenames. Optional, default = 1.") || 
		    !SetParam(argv[1], "ANALYSIS", "ITERATIONS", "Integer. Number of CBA iterations to perform (per time step). Optional, default = 1.") || 
		    !SetParam(argv[1], "ANALYSIS", "COUPLED", "true/false. Whether occupancy is coupled between time steps. Optional, default = false, need OCC_FN when true.") || 
		    !SetParam(argv[1], "ANALYSIS", "DEBUG", "true/false. Whether intermediate files are written each CBA iteration for small datasets (always written for large datasets). Optional, default = false.") || 
		    !SetParam(argv[1], "OUTPUTS", "OCC_FN0", "Filename. Occupancy output filename for time step 0. Additional _FN1..._FNt-1 required for t time steps. Optional if MPC_FN set") || 
		    !SetParam(argv[1], "OUTPUTS", "MPC_FN0", "Filename. MPC output filename for time step 0. Additional _FN1..._FNt-1 required for t time steps. Optional if OCC_FN set") || 
		    !SetParam(argv[1], "OUTPUTS", "EXT_FN0", "Filename. Extinction output filename if using NHA_FN for time step 0. Additional _FN1..._FNt-1 required for t time steps. Optional, default = not written")) 
			return msgErrorStr(-99, "Unable to create parameter file template:", argv[1]);
		return -2;
	}
	testFS.close();
	int returnValue = OccupancyCBA_TS(argv[1]);
	if (returnValue != 0) std::cout << "Error code: " << returnValue << "\n";
	else std::cout << "Complete!" << "\n";
	return returnValue;
}

