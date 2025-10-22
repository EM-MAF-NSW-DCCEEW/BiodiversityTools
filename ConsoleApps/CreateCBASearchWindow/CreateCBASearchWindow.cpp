/*
CreateCBASearchWindow.cpp - Console app for creating a CBA search window
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

#include <fstream>
#include <iostream>

#include "FileUtilsLib.h"
#include "SpatialContextLib.h"

int main(int argc, char* argv[])
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
		if (!SetParam(argv[1], "WINDOW", "LOOKUP_FN", "Filename. Petals look up file name with extension. Optional, Use LOOKUP_FN and TRANSLATE_FN to create petals.") ||
			!SetParam(argv[1], "WINDOW", "TRANSLATE_FN", "Filename. Petals translate file name with extension. Optional, Use LOOKUP_FN and TRANSLATE_FN to create petals..") ||
			!SetParam(argv[1], "WINDOW", "SRCSIZE", "Integer. Source window size. Optional, Use SRCSIZE, DSTSIZE and ZONERATIO to create petals.") ||
			!SetParam(argv[1], "WINDOW", "DSTSIZE", "Integer. Destination window size. Optional, Use SRCSIZE, DSTSIZE and ZONERATIO to create petals.") ||
			!SetParam(argv[1], "WINDOW", "RADIUS", "Floating point number. Radius for window creation. Optional, Use RADIUS, N_RINGS and ZONERATIO to create segments.") ||
			!SetParam(argv[1], "WINDOW", "N_RINGS", "Integer. Number of rings for window creation. Optional, Use RADIUS, N_RINGS and ZONERATIO to create segments.") ||
			!SetParam(argv[1], "WINDOW", "ZONERATIO", "Floating point number. Zone ratio for window scaling. Optional, Use SRCSIZE, DSTSIZE and ZONERATIO to create petals.") ||
			!SetParam(argv[1], "WINDOW", "RADIUSLIST", "Comma-separated list of floating point numbers. Optional, Use RADIUSLIST to create segments.") ||
			!SetParam(argv[1], "WINDOW", "WINDOWGRID_FN", "Filename. Output filename for window grid. Optional, default = not used.") ||
			!SetParam(argv[1], "WINDOW", "PETALDATA_FN", "Filename. Output filename for petal data. Optional, default = not used.")
			)
			return msgErrorStr(-99, "Unable to create parameter file template:", argv[1]);
		return -2;
	}
	testFS.close();
	int returnValue = CreateSearchWindow(argv[1]);
	if (returnValue != 0) std::cout << "Error code: " << returnValue << "\n";
	else std::cout << "Complete!" << "\n";
	return returnValue;
}

