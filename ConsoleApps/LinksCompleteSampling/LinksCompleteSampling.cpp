/*
LinksCompleteSampling.cpp - Console application for Spatial Links with complete sampling 
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

//Notes
//TODO
//Standardise parameter file format

#include <fstream>
#include <iostream>

#include "FileUtilsLib.h"
#include "SpatialLinksLib.h"

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

		if (!SetParam(argv[1], "INPUTS", "CostGrid", "Filename. Cost input grid file name without extension.")) return -1;
		if (!SetParam(argv[1], "INPUTS", "HabGrid", "Filename. Habitat input grid file name without extension. Optional, default = not used.")) return -1;
		if (!SetParam(argv[1], "INPUTS", "SrcHabGrid", "Filename. Source habitat input grid file name without extension. Optional, default = not used.")) return -1;
		if (!SetParam(argv[1], "INPUTS", "DstHabGrid", "Filename. Destination habitat input grid file name without extension. Optional, default = not used.")) return -1;

		if (!SetParam(argv[1], "OUTPUTS", "LinksGrid", "Filename. Links output file name without extension. Optional, default = not created.")) return -1;
		if (!SetParam(argv[1], "OUTPUTS", "CxtGrid", "Filename. Context output file name without extension. Optional, default = not created.")) return -1;

		if (!SetParam(argv[1], "SEARCH", "SearchRadius", "Floating point number. Search radius")) return -1;
		if (!SetParam(argv[1], "SEARCH", "MinHabThreshold", "Floating point number. Minimum habitat threshold")) return -1;
		if (!SetParam(argv[1], "SEARCH", "MaxED", "Floating point number. Maximum effective distance. Optional, default = not used.")) return -1;

		if (!SetParam(argv[1], "DECAY", "I", "Floating point number. Slope of decay function curve")) return -1;
		if (!SetParam(argv[1], "DECAY", "A", "Floating point number. Average movement distance (1/alpha)")) return -1;

		if (!SetParam(argv[1], "OPTIONS", "UseHabWeightLinks", "true/false. Whether to weight links by src/dst habitat. Optional, default = use HabGrid.")) return -1;
		if (!SetParam(argv[1], "OPTIONS", "UseHabWeightCxt", "true/false. Whether to weight context by src habitat. Optional, default = use HabGrid.")) return -1;
		if (!SetParam(argv[1], "OPTIONS", "LinksExponent", "Floating point number. exponent applied to links output")) return -1;
		if (!SetParam(argv[1], "OPTIONS", "nthreads", "Integer. NOT IMPLEMENTED number of threads to use. Optional, default = 1.")) return -1;

		return -2;
	}
	testFS.close();

	int returnValue = 0;
	std::string srcHabFN, dstHabFN;
	if (GetParam(argv[1], "INPUTS", "SrcHabGrid", srcHabFN) && GetParam(argv[1], "INPUTS", "DstHabGrid", dstHabFN)) {
		returnValue = LinksCompleteSamplingAltDstHab(argv[1]);
	}
	else {
		returnValue = LinksCompleteSampling(argv[1]);
	}

	if (returnValue != 0) std::cout << "Error code: " << returnValue << "\n";
	else std::cout << "Complete!" << "\n";
	return returnValue;
}