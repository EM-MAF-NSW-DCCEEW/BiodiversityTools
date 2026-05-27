/*
GDMSimToPixel.cpp - Experimental functions for calculating GDM similarity to a single pixel
Copyright(C) 2023 State of New South Wales and Department of Planning and Environment
Author: Jamie Love, SEI Metrics and Forecasting

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
#include <vector>
//#include "ResilienceDenominator.h"
#include "Petals.h"
#include "FileUtilsLib.h"

//External CUDA functions defined in Common.cuh
extern bool CUDASetDevice(int d);
extern bool CUDAGetGlobalMem(size_t &globMem, int d);

//External CUDA function defined in ResilienceDenominator.cu
extern int CUDAGDMSimToPixel(CBAParams& params);
extern int CUDAGDMSimToPixel(CBAParams& params, uint x, uint y);


///////////////////////////////////////////////////
//Resilience Denominator Sij for samples j
///////////////////////////////////////////////////
int GDMSimToPixel(const std::string & paramFN, CBAParams &params) {
	if (!CUDASetDevice(0)) return msgErrorStr(-10);
	int retval = 0;
	params.maxTxGrids = 32;

	std::string pixelFN, outFN;
	if (!GetParam(paramFN, "INPUTS", "PXL_FN", pixelFN)) return msgErrorStr(-1, "PXL_FN", paramFN);
	if (!GetParam(paramFN, "OUTPUTS", "OUT_FN", outFN)) return msgErrorStr(-1, "OUT_FN", paramFN);

	//Read relisience parameters
	if (!GetParam(paramFN, "RESILIENCE", "NTXGRIDS", params.nTxGrids)) return msgErrorStr(-1, "NTXGRIDS", paramFN);
	if (params.nTxGrids > params.maxTxGrids) return msgErrorStr(-99, "NTXGRIDS exceeds max TX grids");

	if (!ReadGridHeader(pixelFN, params.nCols, params.nRows, params.cellSize, params.noData)) return msgErrorStr(-5, outFN);
	params.nCells = params.nCols * params.nRows;

	params.pixelFS.open(GridFN(pixelFN), std::ios::in | std::ios::binary);
	if (!params.pixelFS.is_open()) return msgErrorStr(-3, outFN);

	params.cxtOutFS.open(GridFN(outFN), std::ios::out | std::ios::binary);
	if (!params.cxtOutFS.is_open()) return msgErrorStr(-4, outFN);
	CopyFS(HeaderFN(pixelFN), HeaderFN(outFN));

	//Open tx data files
	std::string tmpFN;
	for (uint t = 0; t < params.nTxGrids; t++) {
		std::string tStr = toStr(t);

		if (!GetParam(paramFN, "RESILIENCE", "TX" + tStr + "SRC", tmpFN)) return msgErrorStr(-1, "TX" + tStr + "SRC", paramFN);
		if (!CompareGridHeaders(pixelFN, tmpFN)) return msgErrorStr(-7, pixelFN, tmpFN);
		params.txSrcFSs.emplace_back(std::ifstream(GridFN(tmpFN), std::ios::in | std::ios::binary));
		if (!params.txSrcFSs[t].is_open()) return msgErrorStr(-3, tmpFN);

		if (!GetParam(paramFN, "RESILIENCE", "TX" + tStr + "DST", tmpFN)) return msgErrorStr(-1, "TX" + tStr + "DST", paramFN);
		if (!CompareGridHeaders(pixelFN, tmpFN)) return msgErrorStr(-7, pixelFN, tmpFN);
		params.txDstFSs.emplace_back(std::ifstream(GridFN(tmpFN), std::ios::in | std::ios::binary));
		if (!params.txSrcFSs[t].is_open()) return msgErrorStr(-3, tmpFN);
	}

	//Call CUDAResilienceCBA, only bigData
	try {
		retval = CUDAGDMSimToPixel(params);
	}
	catch (std::exception e) {
		return msgErrorStr(-11, e.what());
	}
	return retval;
}
