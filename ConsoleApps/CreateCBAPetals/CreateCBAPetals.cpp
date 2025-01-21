/*
CreateCBAPetals.cpp - Console app for running a context CBA analysis
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

int main(int argc, char *argv[])
{
	int returnValue = -1;
	char fileName[256];
	int srcSize, dstSize;
	float zoneRatio;
	bool createAllFiles = false;
	try {
		if (argc == 2) {
			std::fstream testFS;
			testFS.open(argv[1], std::ios::in);
			if (!testFS.is_open()) {
				std::cout << "\nThe parameter file " << argv[1] << " does not exist and is being created as a template\n\n";
				if (!SetParam(argv[1], "WINDOW", "SRCSIZE", "Integer. Source window size in cells.") || 
				    !SetParam(argv[1], "WINDOW", "DSTSIZE", "Integer. Destination window size in petals.") || 
				    !SetParam(argv[1], "WINDOW", "ZONERATIO", "Floating point number. The area ratio increase of each petal outwards. Must be >= 1.0.") || 
				    !SetParam(argv[1], "WINDOW", "CREATEALLFILES", "true/false. Whether intermediate files are created. Optional, default = false.") || 
				    !SetParam(argv[1], "WINDOW", "PETALS_FN", "FileName. Base output file name with no extension that will be appended to when creating petal files.")) 
					return msgErrorStr(-99, "Unable to create parameter file template:", argv[1]);
				return -2;
			}
			else {
				if (!GetParam(argv[1], "WINDOW", "SRCSIZE", srcSize) || 
				    !GetParam(argv[1], "WINDOW", "DSTSIZE", dstSize) || 
				    !GetParam(argv[1], "WINDOW", "ZONERATIO", zoneRatio) || 
				    !GetParam(argv[1], "WINDOW", "CREATEALLFILES", createAllFiles) || 
				    !GetParam(argv[1], "WINDOW", "PETALS_FN", fileName)) 
					return msgErrorStr(-99, "Unable to read parameters from file:", argv[1]);
			}
		}
		else if (argc < 5) {
			std::cout
				<< "\nUsage: " << argv[0] << " ParameterFile.par \n"
				<< "Or: " << argv[0] << " [FileName] [SrcSize] [DstSize] [ZoneRatio] [CreateAllFiles]\n"
				<< "Or: CreateCBAPetals.exe FileName SrcSize DstSize ZoneRatio CreateAllFiles \n\n"
				<< "If ParameterFile.par exists it is processed, otherwise it is created as a template.\n\n"
				<< "FileName is the base filename with no extension that will be appended to when creating petal files.\n"
				<< "SrcSize and DstSize are the the source window size in cells and destination window size in petals.\n"
				<< "ZoneRatio determines the increase in petal size from each petal zone out to the next.\n"
				<< "CreateAllFiles is true/false determining whether intermediate files are created. Optional, default = false.\n"
				<< "\n"
				<< "SrcSize and DstSize must be odd and DstSize must be less than or equal to SrcSize.\n"
				<< "The maximum allowable DstSize is 31 for a destination window of 31x31 petals.\n"
				<< "ZoneRatio must be greater than or equal to 1.\n"
				<< "The outputs created are FileName_lookup.txt and FileName_translate.txt describing a window configuration used by the CBA.\n"
				<< "If CreateAllFiles equals true or 1, FileName_zoneFactors.txt, FileName_points.asc and FileName_petals.asc are also created.\n"
				<< "FileName_zoneFactors.txt lists the petal zones from the centre out and zone factors used to create petals.\n"
				<< "FileName_points.asc is an ascii grid of petal centroids and FileName_petals.asc is an ascii grid of petals.\n\n";
			return -1;
		}
		else {
			//strcpy_s(fileName, argv[1]);
			strcpy(fileName, argv[1]);
			srcSize = atoi(argv[2]);
			dstSize = atoi(argv[3]);
			zoneRatio = (float)atof(argv[4]);
			if (argc == 6 && (strcmp(argv[5], "true") == 0 || strcmp(argv[5], "1") == 0)) createAllFiles = true;
		}

		returnValue = CreatePetals(fileName, srcSize, dstSize, zoneRatio, createAllFiles);
	}
	catch (std::exception e){
		return msgErrorStr(-11, e.what());
	}

	if (returnValue != 0) std::cout << "Error code: " << returnValue << "\n";
	else std::cout << "Complete!" << "\n";
	return returnValue;
}

