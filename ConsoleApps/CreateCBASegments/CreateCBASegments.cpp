/*
CreateCBASegments.cpp - TODO Description
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


// CreateCBASegments.cpp : Defines the entry point for the console application.
//

#include <fstream>
#include <iostream>
#include <vector>

#include "FileUtilsLib.h"
#include "SpatialContextLib.h"

int main(int argc, char *argv[])
{
	int returnValue = -1;
	//Options
	//Pause after use
	bool useRadList = false;
	char fileName[256];
	int radCells, nRings, nSlices;
	double zoneRatio;
	std::vector<double> radList;

	try {
		if (argc == 2) {
			std::string argv1(argv[1]);
			std::fstream testFS;
			testFS.open(argv[1], std::ios::in);
			if (!testFS.is_open()) {
				std::cout << "\nThe parameter file " << argv[1] << " does not exist and is being created as a template\n\n";
				if (!SetParam(argv1, "WINDOW", "SEGMENTS_FN", "Filename. the base filename used for saving the segments table and grid.") ||
					!SetParam(argv1, "WINDOW", "RAD_CELLS", "Integer. Radius of search area in cells.") ||
				    !SetParam(argv1, "WINDOW", "N_RINGS", "Integer. The number of segment rings.") || 
				    !SetParam(argv1, "WINDOW", "N_SLICES", "Integer. The number of segment slices.") || 
					!SetParam(argv1, "WINDOW", "ZONERATIO", "Floating point number. The area ratio increase of each segment outwards. Must be >= 1.0. Optional, default = use RADLIST.") ||
					!SetParam(argv1, "WINDOW", "RADLIST", "Comma seperated floating point numbers. The radii in cell distance from the focal cell of each ring of segments from smallest to largest. Number of radii must equal N_RINGS.Optional, default = use ZONERATIO."))
					return msgErrorStr(-99, "Unable to create parameter file template:", argv[1]);
				return -2;
			}
			else {
				if (!GetParam(argv[1], "WINDOW", "SEGMENTS_FN", fileName) ||
					!GetParam(argv[1], "WINDOW", "RAD_CELLS", radCells) ||
				    !GetParam(argv[1], "WINDOW", "N_RINGS", nRings) || 
				    !GetParam(argv[1], "WINDOW", "N_SLICES", nSlices)) 
					return msgErrorStr(-99, "Unable to read parameters from file:", argv[1]);

				if (!GetParam(argv[1], "WINDOW", "ZONERATIO", zoneRatio)) {
					if (!GetParam(argv[1], "WINDOW", "RADLIST", radList)) return msgErrorStr(-1, "ZONERATIO or RADLIST", argv[1]);
					useRadList = true;
				}
			}
		}
		else if (argc == 6) {
			strcpy(fileName, argv[1]);
			radCells = atoi(argv[2]);
			nRings = atoi(argv[3]);
			nSlices = atoi(argv[4]);
			zoneRatio = (double)atof(argv[5]);
		} 
		else if (argc > 6) {
			strcpy(fileName, argv[1]);
			radCells = atoi(argv[2]);
			nRings = atoi(argv[3]);
			nSlices = atoi(argv[4]);
			for (int r = 0; r < nRings; r++) radList.push_back(atof(argv[r + 5]));
			useRadList = true;
		}
		else {
			std::cout
				<< "\nUsage: " << argv[0] << " ParameterFile.par \n"
				<< "Or: " << argv[0] << " [Filename] [Radius in cells] [nRings] [nSlices] [Zone ratio]\n"
				<< "Or: " << argv[0] << " [Filename] [Radius in cells] [nRings] [nSlices] [Rad0]...[Radn-1]\n";
			return -1;
		}

		returnValue = useRadList ?
			CreateSegments(radCells, nRings, nSlices, radList, fileName) :
			CreateSegments(radCells, nRings, nSlices, zoneRatio, fileName);
	}

	catch (std::exception e) {
		return msgErrorStr(-11, e.what());
	}

	if (returnValue != 0) std::cout << "Error code: " << returnValue << "\n";
	else std::cout << "Complete!" << "\n";
	return returnValue;
}
