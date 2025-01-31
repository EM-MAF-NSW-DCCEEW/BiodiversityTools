/*
BDIUtils.cpp - TODO Description
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


#include <string>
#include <sstream>
#include <fstream>
#include <cfloat>
#include "FileUtilsLib.h"
#include "BDIUtils.h"

//Get class similarities from sim table if using

bool GetClassSimilarities(const std::string simTabFN, unsigned int nClasses, std::vector<std::vector<double>> & classSim) {
	try {
		if (!FileExists(simTabFN)) return false;
		std::ifstream simTabFS(simTabFN, std::ios::in);
		if (!simTabFS.is_open()) return false;
		classSim.resize(nClasses, std::vector<double>(nClasses, 0.0));
		std::string simLine, simStr;
		for (unsigned int i = 0; i < nClasses; i++) {
			std::getline(simTabFS, simLine);
			std::stringstream simLS(simLine);
			for (unsigned int j = 0; j < nClasses; j++) {
				std::getline(simLS, simStr, ',');
				classSim[i][j] = std::stod(simStr);
			}
		}
		simTabFS.close();
		return true;
	}
	catch (...) { return false; }
}
bool GetClassOHAFromTable(const std::string ohaTabFN, unsigned int nClasses, std::vector<std::vector<double>> & clsOHAInReg) {
	try {
		int cls;
		std::string ohaLine;
		std::ifstream ohaTabFS(ohaTabFN, std::ios::in);
		if (!ohaTabFS.is_open()) return false;
		while (std::getline(ohaTabFS, ohaLine)) {
			if (ohaLine.find(",") > 0 && ohaLine.find(",") < ohaLine.length()) {
				cls = std::stoi(ohaLine.substr(0, ohaLine.find("="))) - 1;
				if (cls >= 0 && cls < signed(nClasses))
					clsOHAInReg[0][cls] = std::stof(ohaLine.substr(ohaLine.find("=") + 1));
			}
		}
		return true;
	}
	catch (...) { return false; }
}

//Finds the smallest power of 10 greater or equal to condMax for scaling condition
bool FindCondDen(float condMax, float &condDen) {
	for (float m = 1.0f; m < 10000000.0f; m *= 10.0f) {
		if (condMax <= m + FLT_EPSILON) {
			condDen = m;
			return true;
		}
	}
	return false;
}