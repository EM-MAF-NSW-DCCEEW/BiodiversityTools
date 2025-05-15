/*
REMP_CBA.cpp - CBA functions for performing REMP analysis
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

#include "Occupancy_CBA.h"
#include "Parameters.h"
#include "Petals.h"
#include "FileUtilsLib.h"

//External CUDA functions defined in Cuda_Common
extern bool CUDASetDevice(int d);
extern bool CUDAGetGlobalMem(size_t &globMem, int d);
extern bool BufferRequired(ull nBands, ull nCells, int d = 0);

//TS REMP
extern int CUDAOccupancyCBA_TS_S(CBAParams &params);
extern int CUDAOccupancyCBA_TS_L(CBAParams &params);

////////////////////////////////////////////////////////////
//Time Series Remp
////////////////////////////////////////////////////////////
int OccupancyCBA_TS(const char * paramFN, CBAParams & params) 
{
	return OccupancyCBA_TS(std::string(paramFN), params);
}

int OccupancyCBA_TS(const std::string & paramFN, CBAParams & p) {
	if (!CUDASetDevice(0)) return msgErrorStr(-10);
	int result = 0;

	p.paramFN = paramFN;
	p.focalPower = 1;

	std::string tmpStr, tmpFN; 
	std::string lookupFN, translateFN, petalsFN;

	//Read required parameters 
	//if (!GetPetalFiles(paramFN, "WINDOW", lookupFN, translateFN, p)) return msgErrorStr(-1, "window parameters", paramFN);
	if (!GetParam(paramFN, "WINDOW", "LOOKUP_FN", lookupFN)) return msgErrorStr(-1, "LOOKUP_FN", paramFN);
	if (!GetParam(paramFN, "WINDOW", "TRANSLATE_FN", translateFN)) return msgErrorStr(-1, "TRANSLATE_FN", paramFN);
	if (!GetParam(paramFN, "SPECIES", "LAMBDA", p.lambda)) return msgErrorStr(-1, "LAMBDA", paramFN);//TODO make optional

	//Read optional parameters
	if (!GetParam(paramFN, "SPECIES", "E", p.e)) p.e = 1.0;
	if (!GetParam(paramFN, "SPECIES", "C", p.c)) p.c = p.e / p.lambda;
	if (!GetParam(paramFN, "ANALYSIS", "MULTFOCAL", p.multFocal)) p.multFocal = false;
	if (!GetParam(paramFN, "ANALYSIS", "FOCALPOWER", p.focalPower)) p.focalPower = 1.0f;
	if (!GetParam(paramFN, "ANALYSIS", "SUMPOWER", p.sumPower)) p.sumPower = 1.0f;
	if (!GetParam(paramFN, "ANALYSIS", "TIMESTEPS", p.rempTimeSteps)) p.rempTimeSteps = 1;
	if (!GetParam(paramFN, "ANALYSIS", "ITERATIONS", p.rempIterations)) p.rempIterations = 1;
	if (!GetParam(paramFN, "ANALYSIS", "MAXPETALPERM", p.useMaxPetalPerm)) p.useMaxPetalPerm = false;
	if (!GetParam(paramFN, "ANALYSIS", "COUPLED", p.rempCoupled)) p.rempCoupled = false;
	if (!GetParam(paramFN, "ANALYSIS", "DEBUG", p.rempWriteEachIteration)) p.rempWriteEachIteration = false;

	//Default to not having optional grids unless they're in the parameter file
	p.haveOccIn = false;
	p.haveMPCIn = false;
	p.haveNHAIn = false;
	p.haveExtIn = false;
	p.haveOccOut = false;
	p.haveMPCOut = false;
	p.haveExtOut = false;

	//OCC and MPC Seeded with optional grid or value
	if (GetParam(paramFN, "INPUTS", "OCC_IN", p.occInFN)) p.haveOccIn = true;
	else if (!GetParam(paramFN, "INPUTS", "OCC_VAL", p.occInitFactor)) p.occInitFactor = 1.0f;

	if (GetParam(paramFN, "INPUTS", "MPC_IN", p.mpcInFN)) p.haveMPCIn = true;
	else if (!GetParam(paramFN, "INPUTS", "MPC_VAL", p.mpcInitFactor)) p.mpcInitFactor = 1.0f;

	//Read all input and output filenames
	for (unsigned int t = 0; t < p.rempTimeSteps; t++) {
		std::string tStr = toStr(t);
		//Inputs
		if (!GetParam(paramFN, "INPUTS", "HAB_FN" + tStr, tmpFN)) return msgErrorStr(-1, "HAB_FN" + tStr, paramFN);
		p.habInFNs.push_back(tmpFN);

		if (!GetParam(paramFN, "INPUTS", "PRM_FN" + tStr, tmpFN)) return msgErrorStr(-1, "PRM_FN" + tStr, paramFN);
		p.prmInFNs.push_back(tmpFN);

		if (GetParam(paramFN, "INPUTS", "NHA_FN" + tStr, tmpFN)) {
			p.haveNHAIn = true;
			p.haveExtIn = false;
			p.nhaInFNs.push_back(tmpFN);
		}

		if (GetParam(paramFN, "INPUTS", "EXT_FN" + tStr, tmpFN) && !p.haveNHAIn) {
			p.haveExtIn = true;
			p.haveNHAIn = false;
			p.extInFNs.push_back(tmpFN);
		}

		//Outputs
		if (GetParam(paramFN, "OUTPUTS", "OCC_FN" + tStr, tmpFN)) {
			p.haveOccOut = true;
			p.occOutFNs.push_back(tmpFN);
		}

		if (GetParam(paramFN, "OUTPUTS", "MPC_FN" + tStr, tmpFN)) {
			p.haveMPCOut = true;
			p.mpcOutFNs.push_back(tmpFN);
		}

		if (GetParam(paramFN, "OUTPUTS", "EXT_FN" + tStr, tmpFN) && p.haveNHAIn) {
			p.haveExtOut = true;
			p.extOutFNs.push_back(tmpFN);
		}

		//Test we have an output and all headers are equal
		if (p.rempCoupled && !p.haveOccOut) return msgErrorStr(-99, "REMP coupled is true and no occupancy output specified");
		if (!p.haveOccOut && !p.haveMPCOut) return msgErrorStr(-99, "No occupancy or MPC output specified");
		if (!ReadGridHeader(p.habInFNs[t], p.nCols, p.nRows, p.cellSize, p.noData)) return msgErrorStr(-5, p.habInFNs[t]);
		if (t > 0 && !CompareGridHeaders(p.habInFNs[0], p.habInFNs[t])) return msgErrorStr(-7, p.habInFNs[0], p.habInFNs[t]);
		if (!CompareGridHeaders(p.habInFNs[t], p.prmInFNs[t])) return msgErrorStr(-7, p.habInFNs[t], p.prmInFNs[t]);
		if (p.haveNHAIn && !CompareGridHeaders(p.habInFNs[t], p.nhaInFNs[t])) return msgErrorStr(-7, p.habInFNs[t], p.nhaInFNs[t]);
		if (p.haveExtIn && !CompareGridHeaders(p.habInFNs[t], p.extInFNs[t])) return msgErrorStr(-7, p.habInFNs[t], p.extInFNs[t]);
		if (p.haveOccOut && !CopyFS(HeaderFN(p.habInFNs[t]), HeaderFN(p.occOutFNs[t]), true)) return msgErrorStr(-6, HeaderFN(p.occOutFNs[t]));
		if (p.haveMPCOut && !CopyFS(HeaderFN(p.habInFNs[t]), HeaderFN(p.mpcOutFNs[t]), true)) return msgErrorStr(-6, HeaderFN(p.mpcOutFNs[t]));
	}

	//Create petal table
	if (!CreatePetalData(translateFN, lookupFN, RemoveFileExt(paramFN), p)) return msgErrorStr(-9, petalsFN);

	//Set nCells
	p.nCells = p.nCols * p.nRows;
	
	//Call CUDA_REMP_TS_CBA
	try {
		//result = BufferRequired(6ULL, p.nCells) ? CUDAOccupancyCBA_TS_L(p) : CUDAOccupancyCBA_TS_S(p);

		//result = CUDAOccupancyCBA_TS_L(p);

		result = CUDAOccupancyCBA_TS_S(p);


	}
	catch (std::exception e) {
		return msgErrorStr(-11, e.what());
	}
	return result;
}

