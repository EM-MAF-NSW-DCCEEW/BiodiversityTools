/*
Benefits_CBA.cpp - CBA functions for performing BFT benefits analysis
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

#include "Benefits_CBA.h"
#include "Petals.h"
#include "FileUtilsLib.h"

//External CUDA functions defined in Common.cuh
extern bool CUDASetDevice(int d);
extern bool CUDAGetGlobalMem(size_t &globMem, int d);
extern bool BufferRequired(ull nBands, ull nCells, int d = 0);

//External CUDA functions defined in Benefits_CBA.cu
extern int CUDABenefitCBA_S(CBAParams &params);
extern int CUDABenefitCBA_L(CBAParams &params);

///////////////////////////////////////////////////
//Benefit CBA
///////////////////////////////////////////////////
int BenefitsCBA(const char *paramFN, CBAParams &params) 
{ 
	return BenefitsCBA(std::string(paramFN), params);
}

int BenefitsCBA(const std::string & paramFN, CBAParams &params)
{
	if (!CUDASetDevice(0)) return msgErrorStr(-10);
	int result = 0;

	//FileNames
	std::string habFN, prmFN, mbvFN, altHabFN, altPrmFN, benFN;
	std::string pLookupFN, pTranslateFN, petalsFN;

	//Read required parameters
	if (!GetParam(paramFN, "INPUTS", "HAB_FN", habFN)) return msgErrorStr(-1, "HAB_FN", paramFN);
	if (!GetParam(paramFN, "INPUTS", "PRM_FN", prmFN)) return msgErrorStr(-1, "PRM_FN", paramFN);
	if (!GetParam(paramFN, "INPUTS", "MBV_FN", mbvFN)) return msgErrorStr(-1, "MBV_FN", paramFN);
	if (!GetParam(paramFN, "OUTPUTS", "BEN_FN", benFN)) return msgErrorStr(-1, "BEN_FN", paramFN);
	if (!GetParam(paramFN, "WINDOW", "LOOKUP_FN", pLookupFN)) return msgErrorStr(-1, "LOOKUP_FN", paramFN);
	if (!GetParam(paramFN, "WINDOW", "TRANSLATE_FN", pTranslateFN)) return msgErrorStr(-1, "TRANSLATE_FN", paramFN);

	//Read optional parameters
	if (!GetParam(paramFN, "ANALYSIS", "MULTFOCAL", params.multFocal)) params.multFocal = true;
	if (!GetParam(paramFN, "ANALYSIS", "FOCALPOWER", params.focalPower)) params.focalPower = 1.0f;
	if (!GetParam(paramFN, "ANALYSIS", "SUMPOWER", params.sumPower)) params.sumPower = 1.0f;

	//Read required alt hab and prm parameters
	if (!GetParam(paramFN, "INPUTS", "ALTGRIDS", params.useAltGrids)) return msgErrorStr(-1, "ALTGRIDS", paramFN);
	if (params.useAltGrids) {
		if (!GetParam(paramFN, "INPUTS", "ALTHAB", altHabFN)) return msgErrorStr(-1, "ALTHAB", paramFN);
		if (!GetParam(paramFN, "INPUTS", "ALTPRM", altPrmFN)) return msgErrorStr(-1, "ALTPRM", paramFN);
	}
	else {
		if (!GetParam(paramFN, "INPUTS", "ALTHAB", params.altHab)) return msgErrorStr(-1, "ALTHAB", paramFN);
		if (!GetParam(paramFN, "INPUTS", "ALTPRM", params.altPrm)) return msgErrorStr(-1, "ALTPRM", paramFN);
	}

	//TODO Test this to remove ALTGRIDS parameter
	//if (!(GetParam(paramFN, "INPUTS", "ALTHAB", params.altHab) || 
	//	(GetParam(paramFN, "INPUTS", "ALTHAB", altHabFN) && FileExists(GridFN(altHabFN))))) 
	//	return msgErrorStr(-1, "ALTHAB", paramFN);
	//if (!(GetParam(paramFN, "INPUTS", "ALTPRM", params.altPrm) || 
	//	(GetParam(paramFN, "INPUTS", "ALTPRM", altPrmFN) && FileExists(GridFN(altPrmFN))))) 
	//	return msgErrorStr(-1, "ALTPRM", paramFN);


	//Open data files
	params.habInFS.open(GridFN(habFN), std::ios::in | std::ios::binary);
	params.prmInFS.open(GridFN(prmFN), std::ios::in | std::ios::binary);
	params.mbvInFS.open(GridFN(mbvFN), std::ios::in | std::ios::binary);
	params.cxtOutFS.open(GridFN(benFN), std::ios::out | std::ios::binary);
	if (!params.habInFS.is_open()) return msgErrorStr(-3, habFN);
	if (!params.prmInFS.is_open()) return msgErrorStr(-3, prmFN);
	if (!params.mbvInFS.is_open()) return msgErrorStr(-3, mbvFN);
	if (!params.cxtOutFS.is_open()) return msgErrorStr(-4, benFN);

	if (params.useAltGrids) {
		params.altHabInFN.open(GridFN(altHabFN), std::ios::in | std::ios::binary);
		params.altPrmInFS.open(GridFN(altPrmFN), std::ios::in | std::ios::binary);
		if (!params.altHabInFN.is_open()) return msgErrorStr(-3, altHabFN);
		if (!params.altPrmInFS.is_open()) return msgErrorStr(-3, altPrmFN);
	}

	//Read and test header files
	if (!ReadGridHeader(habFN, params.nCols, params.nRows, params.cellSize, params.noData)) return msgErrorStr(-5, habFN);
	if (!CompareGridHeaders(habFN, prmFN)) return msgErrorStr(-7, habFN, prmFN);
	if (!CompareGridHeaders(habFN, mbvFN)) return msgErrorStr(-7, habFN, mbvFN);

	//Test altHab header values
	if (params.useAltGrids) {
		if (!CompareGridHeaders(habFN, altHabFN)) return msgErrorStr(-7, habFN, altHabFN);
		if (!CompareGridHeaders(habFN, altPrmFN)) return msgErrorStr(-7, habFN, altPrmFN);
	}

	//Create output header file
	if (!CopyFS(HeaderFN(habFN), HeaderFN(benFN), true)) return msgErrorStr(-8, habFN, benFN);

	//Create petal table
	if (!CreatePetalData(pTranslateFN, pLookupFN, benFN, params)) return msgErrorStr(-9, petalsFN);

	//Set nCells
	params.nCells = params.nCols * params.nRows;

	//Call Benefit CBA
	try {
		result = BufferRequired(6ULL, params.nCells) ? CUDABenefitCBA_L(params) : CUDABenefitCBA_S(params);
	}
	catch (std::exception e) {
		return msgErrorStr(-11, e.what());
	}
	return result;
}

