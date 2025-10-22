/*
ContinuousMBV.cpp - Experimental functions for calculating marginal biodiversity values from continuous GDM data
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
#include "ContinuousMBV.h"
#include "Petals.h"
#include "FileUtilsLib.h"

#include <climits>
#include <cmath>

//External CUDA functions defined in Common.cuh
extern bool CUDASetDevice(int d);
extern bool CUDAGetGlobalMem(size_t &globMem, int d);

//External CUDA function defined in continuousBDI.cu
extern int CUDAContinuousMBV(CBAParams &params);

///////////////////////////////////////////////////
//Resilience CBA
///////////////////////////////////////////////////
int ContinuousMBV(const std::string & paramFN, CBAParams &params) {
	if (!CUDASetDevice(0)) return msgErrorStr(-10);
	int retval = 0;

	std::string habFN;
	std::string outFN;
	if (!GetParam(paramFN, "RESILIENCE", "HAB_FN", habFN)) return msgErrorStr(-1, "HAB_FN", paramFN);
	if (!GetParam(paramFN, "OUTPUTS", "OUT_FN", outFN)) return msgErrorStr(-1, "OUT_FN", paramFN);


	//Read relisience parameters
	params.sampleRate = 0;
	params.maxSamples = 0;
	GetParam(paramFN, "RESILIENCE", "SAMPLERATE", params.sampleRate);
	GetParam(paramFN, "RESILIENCE", "MAXSAMPLES", params.maxSamples);
	if (!GetParam(paramFN, "RESILIENCE", "NTXGRIDS", params.nTxGrids)) return msgErrorStr(-1, "NTXGRIDS", paramFN);

	params.maxTxGrids = 32;

	//Open data files
	params.habInFS.open(GridFN(habFN), std::ios::in | std::ios::binary);
	if (!params.habInFS.is_open()) return msgErrorStr(-3, habFN);

	params.numOutFS.open(GridFN(outFN + "_num"), std::ios::out | std::ios::binary);
	if (!params.numOutFS.is_open()) return msgErrorStr(-4, outFN + "_num");

	params.denOutFS.open(GridFN(outFN + "_den"), std::ios::out | std::ios::binary);
	if (!params.denOutFS.is_open()) return msgErrorStr(-4, outFN + "_den");

	params.txDstSamplesFS.open(outFN + ".samples", std::ios::out);
	if (!params.txDstSamplesFS.is_open()) return msgErrorStr(-4, outFN + ".samples");

	params.sampleGridFS.open(GridFN(outFN + "_samples"), std::ios::out | std::ios::binary);
	if (!params.txDstSamplesFS.is_open()) return msgErrorStr(-4, GridFN(outFN + "_samples"));


	//Count number of data cells in hab grid
	if (!ReadGridHeader(habFN, params.nCols, params.nRows, params.cellSize, params.noData)) return msgErrorStr(-5, outFN);
	params.nCells = params.nCols * params.nRows;

	ull nDataCells = 0;

	//Open tx data files
	std::string tmpFN;
	for (uint t = 0; t < params.nTxGrids; t++) {
		std::string tStr = toStr(t);

		if (!GetParam(paramFN, "RESILIENCE", "TX" + tStr + "SRC", tmpFN)) return msgErrorStr(-1, "TX" + tStr + "SRC", paramFN);
		//if (t > 0 && !CompareGridHeaders(habFN, tmpFN)) return msgErrorStr(-7, habFN, tmpFN);
		//params.txSrcFSs.emplace_back(std::ifstream(GridFN(tmpFN), std::ios::in | std::ios::binary));
		//if (!params.txSrcFSs[t].is_open()) return msgErrorStr(-3, tmpFN);
		params.txSrcGSs.emplace_back(igstream(tmpFN));
		if (!params.txSrcGSs[t].is_open()) return msgErrorStr(-3, tmpFN);

		if (!GetParam(paramFN, "RESILIENCE", "TX" + tStr + "DST", tmpFN)) return msgErrorStr(-1, "TX" + tStr + "DST", paramFN);
		//if (t > 0 && !CompareGridHeaders(habFN, tmpFN)) return msgErrorStr(-7, habFN, tmpFN);
		//params.txDstFSs.emplace_back(std::ifstream(GridFN(tmpFN), std::ios::in | std::ios::binary));
		//if (!params.txDstFSs[t].is_open()) return msgErrorStr(-3, tmpFN);
		params.txDstGSs.emplace_back(igstream(tmpFN));
		if (!params.txDstGSs[t].is_open()) return msgErrorStr(-3, tmpFN);

		////Create output header file
		//if (t == 0) {
		//	if (!CopyFS(HeaderFN(tmpFN), HeaderFN(outFN + "_num"), true)) return msgErrorStr(-8, HeaderFN(tmpFN), HeaderFN(outFN + "_num"));
		//	if (!CopyFS(HeaderFN(tmpFN), HeaderFN(outFN + "_den"), true)) return msgErrorStr(-8, HeaderFN(tmpFN), HeaderFN(outFN + "_den"));
		//	if (!CopyFS(HeaderFN(tmpFN), HeaderFN(outFN + "_samples"), true)) return msgErrorStr(-8, HeaderFN(tmpFN), HeaderFN(outFN + "_samples"));
		//	CountDataCells(GridFN(tmpFN), nDataCells);
		//}
	}

	//Create output header file
	if (!CopyFS(HeaderFN(habFN), HeaderFN(outFN + "_num"), true)) return msgErrorStr(-8, HeaderFN(tmpFN), HeaderFN(outFN + "_num"));
	if (!CopyFS(HeaderFN(habFN), HeaderFN(outFN + "_den"), true)) return msgErrorStr(-8, HeaderFN(tmpFN), HeaderFN(outFN + "_den"));
	if (!CopyFS(HeaderFN(habFN), HeaderFN(outFN + "_samples"), true)) return msgErrorStr(-8, HeaderFN(tmpFN), HeaderFN(outFN + "_samples"));
	
	CountDataCells(GridFN(habFN), nDataCells);
	if (nDataCells > ull(UINT_MAX)) return msgError(-99, "nDataCells overflows uint. Grid sizes over 2^32 not supported");
	
	if (params.maxSamples > 0) params.sampleRate = uint(nDataCells) / params.maxSamples;
	else if (params.sampleRate > 0) params.maxSamples = uint(nDataCells) / params.sampleRate;
	else return -1;
	params.sampleStep = int(sqrt(params.sampleRate));
	std::cout << "Number of cells: \t\t" << params.nCells << "\n";
	std::cout << "Data cells: \t\t" << nDataCells << "\n";
	std::cout << "Sample Rate: \t\t" << params.sampleRate << "\n";
	std::cout << "Max Samples: \t\t" << params.maxSamples << "\n";
	std::cout << "Sample Step: \t\t" << params.sampleStep << "\n";

	//Call CUDAResilienceCBA, only bigData
	try {
		retval = CUDAContinuousMBV(params);
	}
	catch (std::exception e) {
		return msgErrorStr(-11, e.what());
	}
	return retval;
}
