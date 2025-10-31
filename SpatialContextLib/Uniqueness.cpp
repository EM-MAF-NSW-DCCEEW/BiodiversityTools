/*
Uniqueness.cpp - Functions for calculating uniqueness from continuous GDM data
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
#include "uniqueness.h"
#include "Petals.h"
#include "FileUtilsLib.h"

#include <climits>
#include <cmath>

//External CUDA functions defined in Common.cuh
extern bool CUDASetDevice(int d);
extern bool CUDAGetGlobalMem(size_t& globMem, int d);

//External CUDA function defined in ResilienceDenominator.cu
extern int CUDAUniqueness(CBAParams& params);

//Get the number of data cells in a flt grid
bool GetnDataCells(const std::string& FN, ull& nDataCells)
{
	try {
		uint nCols, nRows;
		double cellSize;
		float noData;
		ull count = 0;
		igstream gridGS(FN);
		if (!gridGS.is_open()) return false;
		if (!gridGS.getHeader(nCols, nRows, cellSize, noData)) return false;

		std::vector<float> data(nCols);
		
		for (uint r = 0; r < nRows; r++) {
			gridGS.read((char*)data.data(), nCols * sizeof(float));
			for (uint c = 0; c < nCols; c++) {
				if (data[c] == noData) continue;
				if (ULLONG_MAX - count < 1ULL) return false;
				++count;
			}
		}
		nDataCells = count;
		return true;
	}
	catch (...) { return false; }
}


///////////////////////////////////////////////////
//Uniqueness ~Sij for samples j
///////////////////////////////////////////////////
int Uniqueness(const std::string& paramFN, CBAParams& params) {
	if (!CUDASetDevice(0)) return msgErrorStr(-10);
	int retval = 0;

	//Get output file names
	std::string outputGridFN;
	if (!GetParam(paramFN, "OUTPUTS", "OUT_FN", outputGridFN)) return msgErrorStr(-1, "OUT_FN", paramFN);
	std::string txDstSamplesFN = RemoveFileExt(outputGridFN) + "_samples.csv";
	std::string sampleGridFN = RemoveFileExt(outputGridFN) + "_samples.tif";

	//Read relisience parameters
	params.sampleRate = 0;
	params.maxSamples = 0;
	GetParam(paramFN, "RESILIENCE", "SAMPLERATE", params.sampleRate);
	GetParam(paramFN, "RESILIENCE", "MAXSAMPLES", params.maxSamples);
	if (!GetParam(paramFN, "RESILIENCE", "NTXGRIDS", params.nTxGrids)) return msgErrorStr(-1, "NTXGRIDS", paramFN);

	params.maxTxGrids = 32;

	//Open tx data files
	std::string tmpFN;
	for (uint t = 0; t < params.nTxGrids; t++) {
		std::string tStr = toStr(t);

		if (!GetParam(paramFN, "RESILIENCE", "TX" + tStr + "SRC", tmpFN)) return msgErrorStr(-1, "TX" + tStr + "SRC", paramFN);
		params.txSrcGSs.emplace_back(igstream(tmpFN));
		if (!params.txSrcGSs[t].is_open()) return msgErrorStr(-3, tmpFN);
		if (t > 0 && !(params.txSrcGSs[t].compareHeader(params.txSrcGSs[0]))) return msgErrorStr(-7, params.txSrcGSs[t].fileName(), params.txSrcGSs[0].fileName());

		if (!GetParam(paramFN, "RESILIENCE", "TX" + tStr + "DST", tmpFN)) return msgErrorStr(-1, "TX" + tStr + "DST", paramFN);
		params.txDstGSs.emplace_back(igstream(tmpFN));
		if (!params.txDstGSs[t].is_open()) return msgErrorStr(-3, tmpFN);
		if (t > 0 && !(params.txDstGSs[t].compareHeader(params.txDstGSs[0]))) return msgErrorStr(-7, params.txDstGSs[t].fileName(), params.txDstGSs[0].fileName());
	}

	//Get header info from first txDst grid and set output headers
	if (!params.txDstGSs[0].getHeader(params.nCols, params.nRows, params.cellSize, params.noData)) return msgErrorStr(-5, params.txDstGSs[0].fileName());;
	
	//Open output files
	params.outUniquenessGS.open(outputGridFN);
	params.outUniquenessGS.copyHeader(params.txDstGSs[0]);
	params.outUniquenessGS.copyProjection(params.txDstGSs[0]);
	if (!params.outUniquenessGS.is_open()) return msgErrorStr(-4, outputGridFN);

	params.sampleGridGS.open(sampleGridFN);
	params.sampleGridGS.copyHeader(params.txDstGSs[0]);
	params.sampleGridGS.copyProjection(params.txDstGSs[0]);
	if (!params.sampleGridGS.is_open()) return msgErrorStr(-4, sampleGridFN);

	params.txDstSamplesFS.open(txDstSamplesFN, std::ios::out);
	if (!params.txDstSamplesFS.is_open()) return msgErrorStr(-4, txDstSamplesFN);

	//Calculate number of cells
	params.nCells = params.nCols * params.nRows;

	//Count number of data cells in first txDst grid
	ull nDataCells = 0;
	GetnDataCells(params.txDstGSs[0].fileName(), nDataCells);
	if (nDataCells == 0) return msgError(-99, "No data cells found in TX DST grid");
	if (nDataCells > ull(UINT_MAX)) return msgError(-99, "nDataCells overflows uint. Grid sizes over 2^32 not supported");
	
	//Determine sample rate/ max samples
	if (params.maxSamples > 0) params.sampleRate = uint(nDataCells) / params.maxSamples;
	else if (params.sampleRate > 0) params.maxSamples = uint(nDataCells) / params.sampleRate;
	else return msgError(-99, "Neither SAMPLERATE or MAXSAMPLES greater than 0");

	params.sampleStep = int(std::sqrt(params.sampleRate));
	std::cout << "Number of cells: \t\t" << params.nCells << "\n";
	std::cout << "Number of data cells: \t\t" << nDataCells << "\n";
	std::cout << "Sample Rate: \t\t" << params.sampleRate << "\n";
	std::cout << "Max Samples: \t\t" << params.maxSamples << "\n";
	std::cout << "Sample Step: \t\t" << params.sampleStep << "\n";

	//Call CUDA Uniqueness function
	try {
		retval = CUDAUniqueness(params);
	}
	catch (std::exception e) {
		return msgErrorStr(-11, e.what());
	}
	return retval;
}
