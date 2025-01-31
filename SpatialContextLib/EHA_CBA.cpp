/*
EHA_CBA.cpp - CBA functions for performing Effective Habitat Area analysis
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
#include "EHA_CBA.h"
#include "Petals.h"
#include "FileUtilsLib.h"

//External CUDA functions defined in Common.cuh
extern bool CUDASetDevice(int d);
extern bool CUDAGetGlobalMem(size_t &globMem, int d);
extern bool BufferRequired(ull nBands, ull nCells, int d = 0);

//External CUDA functions defined in Context_CBA.cu
extern int CUDA_EHA_CBA_S(CBAParams &params);
extern int CUDA_EHA_CBA_L(CBAParams &params);

///////////////////////////////////////////////////
//EHA CBA
///////////////////////////////////////////////////
int EHA_CBA(const char *paramFN, CBAParams &params) 
{
	return EHA_CBA(std::string(paramFN), params);
}

int EHA_CBA(const std::string & paramFN, CBAParams &params) {

	if (!CUDASetDevice(0)) return msgErrorStr(-10);
	int result = 0;

	//FileNames
	std::string habFN, prmFN, cxtFN;
	std::string pLookupFN, pTranslateFN, petalsFN;

	//Read required parameters
	if (!GetParam(paramFN, "INPUTS", "HAB_FN", habFN)) return msgErrorStr(-1, "HAB_FN", paramFN);
	if (!GetParam(paramFN, "INPUTS", "PRM_FN", prmFN)) return msgErrorStr(-1, "PRM_FN", paramFN);
	if (!GetParam(paramFN, "OUTPUTS", "EHA_FN", cxtFN)) return msgErrorStr(-1, "EHA_FN", paramFN);
	if (!GetParam(paramFN, "WINDOW", "LOOKUP_FN", pLookupFN)) return msgErrorStr(-1, "LOOKUP_FN", paramFN);
	if (!GetParam(paramFN, "WINDOW", "TRANSLATE_FN", pTranslateFN)) return msgErrorStr(-1, "TRANSLATE_FN", paramFN);

	//Read optional parameters
	if (!GetParam(paramFN, "ANALYSIS", "MULTFOCAL", params.multFocal)) params.multFocal = false;
	if (!GetParam(paramFN, "ANALYSIS", "FOCALPOWER", params.focalPower)) params.focalPower = 1.0f;
	if (!GetParam(paramFN, "ANALYSIS", "SUMPOWER", params.sumPower)) params.sumPower = 1.0f;

	//Get max Cond and Perm
	GetGridMaxPow10(habFN, params.maxCond);
	GetGridMax(prmFN, params.maxPerm);

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

	//Create output header file
	if (!CopyFS(HeaderFN(habFN), HeaderFN(cxtFN), true)) return msgErrorStr(-8, habFN, cxtFN);

	//Set nCells
	params.nCells = params.nCols * params.nRows;

	//Create petal table
	if (!CreatePetalData(pTranslateFN, pLookupFN, cxtFN, params)) return msgErrorStr(-9, petalsFN);

	//Call Factor or Context CBA
	try {
		result = BufferRequired(3ULL, params.nCells) ? CUDA_EHA_CBA_L(params) : CUDA_EHA_CBA_S(params);
	}
	catch (std::exception e) {
		return msgErrorStr(-11, e.what());
	}
	return result;
}
