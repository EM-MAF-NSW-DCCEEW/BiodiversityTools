/*
Context_CBA.cpp - CBA functions for performing BFT context analysis
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

#include <stdio.h>
#include <fstream>
#include <iostream>
#include "Context_CBA.h"
#include "Petals.h"
#include "FileUtilsLib.h"

//TODO update to use window parameters or segments rather than lookup and translate files


//External CUDA functions defined in Common.cuh
extern bool CUDASetDevice(int d);
extern bool CUDAGetGlobalMem(size_t &globMem, int d);
extern bool BufferRequired(ull nBands, ull nCells, int d = 0);

//External CUDA functions defined in Context_CBA.cu
extern int CUDAContextCBA_S(CBAParams &params);
extern int CUDAContextCBA_L(CBAParams &params);

//extern int CUDAContextProbCBA_S(CBAParams &params);
//extern int CUDAContextProbCBA_L(CBAParams &params);

//External CUDA functions defined in Factor_CBA.cu
extern int CUDAFactorCBA_S(CBAParams &params);
extern int CUDAFactorCBA_L(CBAParams &params);

///////////////////////////////////////////////////
//Context CBA
///////////////////////////////////////////////////
int ContextCBA(const char *paramFN, CBAParams &params) 
{
	return ContextCBA(std::string(paramFN), params);
}

int ContextCBA(const std::string & paramFN, CBAParams &params)
{
	if (!CUDASetDevice(0)) return msgErrorStr(-10);
	int result = 0;

	//FileNames
	std::string habFN, prmFN, facFN, cxtFN;
	std::string pLookupFN, pTranslateFN, petalsFN;

	//Resilience defaults to adaptive buffering
	bool useSegments = false;
	bool useRadList = false;

	int radCells, nRings, nSlices;
	double zRatio;
	std::vector<double> radList;

	//Read required parameters
	if (!GetParam(paramFN, "INPUTS", "HAB_FN", habFN)) return msgErrorStr(-1, "HAB_FN", paramFN);
	if (!GetParam(paramFN, "INPUTS", "PRM_FN", prmFN)) return msgErrorStr(-1, "PRM_FN", paramFN);
	if (!GetParam(paramFN, "OUTPUTS", "NHA_FN", cxtFN)) return msgErrorStr(-1, "NHA_FN", paramFN);
	//if (!GetParam(paramFN, "WINDOW", "LOOKUP_FN", pLookupFN)) return msgErrorStr(-1, "LOOKUP_FN", paramFN);
	//if (!GetParam(paramFN, "WINDOW", "TRANSLATE_FN", pTranslateFN)) return msgErrorStr(-1, "TRANSLATE_FN", paramFN);

	//Read window parameters
	if (GetParam(paramFN, "WINDOW", "LOOKUP_FN", pLookupFN) &&
		GetParam(paramFN, "WINDOW", "TRANSLATE_FN", pTranslateFN))
		useSegments = false;
	else if (GetParam(paramFN, "WINDOW", "RAD_CELLS", radCells) &&
		GetParam(paramFN, "WINDOW", "N_RINGS", nRings) &&
		GetParam(paramFN, "WINDOW", "N_SLICES", nSlices)) {
		useSegments = true;
		if (!GetParam(paramFN, "WINDOW", "ZONERATIO", zRatio)) {
			if (!GetParam(paramFN, "WINDOW", "RADLIST", radList)) return msgErrorStr(-1, "ZONERATIO or RADLIST", paramFN);
			useRadList = true;
		}
	}
	else return msgErrorStr(-1, "window parameters", paramFN);





	//Read optional parameters
	params.haveFacIn = GetParam(paramFN, "INPUTS", "FAC_FN", facFN);
	if (!GetParam(paramFN, "ANALYSIS", "MULTFOCAL", params.multFocal)) params.multFocal = false;
	if (!GetParam(paramFN, "ANALYSIS", "FOCALPOWER", params.focalPower)) params.focalPower = 1.0f;
	if (!GetParam(paramFN, "ANALYSIS", "SUMPOWER", params.sumPower)) params.sumPower = 1.0f;

	if (!GetParam(paramFN, "DISPERSAL", "AVGHAB", params.doAveragePetalHabitat)) params.doAveragePetalHabitat = false;
	if (!GetParam(paramFN, "DISPERSAL", "DISPPROB", params.doDispersalProbability)) params.doDispersalProbability = false;
	if (!GetParam(paramFN, "DISPERSAL", "DISPYEARS", params.dispersalYears)) params.dispersalYears= 50.0f;



	//Open data files
	params.habInFS.open(GridFN(habFN), std::ios::in | std::ios::binary);
	params.prmInFS.open(GridFN(prmFN), std::ios::in | std::ios::binary);
	params.cxtOutFS.open(GridFN(cxtFN), std::ios::out | std::ios::binary);
	if (!params.habInFS.is_open()) return msgErrorStr(-3, habFN);
	if (!params.prmInFS.is_open()) return msgErrorStr(-3, prmFN);
	if (!params.cxtOutFS.is_open()) return msgErrorStr(-4, cxtFN);

	//Read and test header files
	if (!ReadGridHeader(habFN, params.nCols, params.nRows, params.cellSize, params.noData)) return msgErrorStr(-5, habFN);
	if (!CompareGridHeaders(habFN, prmFN)) return msgErrorStr(-7, habFN, prmFN);

	//Use factor grid
	if (params.haveFacIn) {
		params.facInFS.open(GridFN(facFN), std::ios::in | std::ios::binary);
		if (!params.facInFS.is_open()) return msgErrorStr(-3, facFN);
		if (!CompareGridHeaders(habFN, facFN)) return msgErrorStr(-7, habFN, facFN);
	}

	//Create output header file
	if (!CopyFS(HeaderFN(habFN), HeaderFN(cxtFN), true)) return msgErrorStr(-8, habFN, cxtFN);

	//Set nCells
	params.nCells = params.nCols * params.nRows;

	//Create petal table
	//if (!CreatePetalData(pTranslateFN, pLookupFN, cxtFN, params)) return msgErrorStr(-9, petalsFN);

	//Create petal table
	if (useSegments) {
		int retval = useRadList ?
			CreateSegments(radCells, nRings, nSlices, radList, petalsFN, params) :
			CreateSegments(radCells, nRings, nSlices, zRatio, petalsFN, params);
		if (retval != 0) return msgErrorStr(-99, "Unable to create segments:", petalsFN);
	}
	else {
		if (!CreatePetalData(pTranslateFN, pLookupFN, cxtFN, params)) return msgErrorStr(-9, petalsFN);
	}

	//Call CUDAFactorCBA or CUDAContextCBA
	try {
		if (params.haveFacIn)
			result = BufferRequired(3ULL, params.nCells) ? CUDAFactorCBA_L(params) : CUDAFactorCBA_S(params);
		else
			result = BufferRequired(5ULL, params.nCells) ? CUDAContextCBA_L(params) : CUDAContextCBA_S(params);
	}
	catch (std::exception e) {
		return msgErrorStr(-11, e.what());
	}
	return result;
}

