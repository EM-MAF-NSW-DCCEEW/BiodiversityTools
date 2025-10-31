/*
Parameters.h defines CBAParams class used for storing CBA parameters,filenames etc
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

#pragma once
#ifndef _CUDA_CBA_PARAMETERS_H_
#define _CUDA_CBA_PARAMETERS_H_
//#include <cuda_runtime.h>
//#include <device_launch_parameters.h
#include "vector_types.h"
#include <vector>
#include <map>
#include <list>
#include <fstream>
#include "FileUtilsLib.h"
#include <memory>

using uint = unsigned int;
using ull = unsigned long long;

#ifdef _WIN32
#   define STDCALL __stdcall
#else //LINUX
#   define STDCALL 
#endif

//CUDA CBA Parameters contains Filestreams so no copy
struct CBAParams {
	
	//Callback function for BFT
	bool(STDCALL *fptr)(char *p, int i);

	//Parameter File Name
	std::string paramFN;

	//Grid parameters
	uint nCols, nRows;
	ull nCells;
	float noData; 
	double cellSize;

	//Window parameters
	uint petalRows, petalCols, nPetals, fOffset, firstReduction;
	std::unique_ptr<int2[]> petalPtr;

	//CBA Options
	float sumPower, focalPower;
	bool multFocal;

	//Context CBA
	bool haveFacIn;
	std::ifstream habInFS;
	std::ifstream prmInFS;
	std::ifstream facInFS;
	std::ofstream cxtOutFS;
	
	//Benefits CBA
	float altHab, altPrm;
	bool useAltGrids;
	std::ifstream mbvInFS;
	std::ifstream altHabInFN;
	std::ifstream altPrmInFS;

	//EHA CBA
	float maxCond, maxPerm;
	
	//REMP CBA
	uint rempTimeSteps = 0, rempIterations = 0;
	float occInitFactor, mpcInitFactor, mpcDen;
	double c, e, lambda;
	bool haveOccIn, haveMPCIn, haveNHAIn, haveExtIn, haveExtOut, haveOccOut, haveMPCOut, rempWriteEachIteration, rempCoupled, useMaxPetalPerm;
	std::ifstream nhaInFS;
	std::ifstream extInFS;
	std::ifstream occInFS;
	std::ifstream mpcInFS;
	std::ofstream occPuOutFS;
	std::ofstream occSiOutFS;
	std::ofstream mpcOutFS;
	std::ofstream extOutFS;
	std::fstream mpcTmpFS;
	std::string occInFN, mpcInFN, occPuOutFN, occSiOutFN, mpcOutFN, tmpFN;
	std::vector<std::string> habInFNs, prmInFNs, nhaInFNs, extInFNs, occOutFNs, mpcOutFNs, extOutFNs;

	//RESILIENCE CBA
	std::ofstream numOutFS, denOutFS;
	uint nTxGrids = 0;
	uint colStep = 32;
	std::vector<std::ifstream> txSrcFSs, txDstFSs;
	//Updated to use gstreams and read tif transform grids - Already declared below
	//std::vector<igstream> txSrcGSs, txDstGSs;

	//RESILIENCE MULTI_DISPERSAL
	uint nDispersals;
	std::vector<float2> perms;
	std::vector<std::ofstream> outFSs;

	//RESILIENCE DENOMINATOR & CONTINUOUS MBV
	uint maxSamples;
	uint sampleRate;
	uint sampleStep;
	uint maxTxGrids;
	float sampleRateF;
	std::ofstream txDstSamplesFS;
	std::ofstream sampleGridFS;

	//RESILIENCE TESTING DISPERSAL PROBABILITY
	bool doAveragePetalHabitat;
	bool doDispersalProbability;
	float dispersalYears;

	//RESILIENCE TESTING GAUSSIAN DISPERSAL
	uint nGenerations;
	float habFactor;
	float gaussConst;
	std::vector<float> lambdas;

	//GDM SIMILARITY TO CLASSES OR SPECIES OBSERVATIONS
	std::ifstream clsFS, obsFS;
	std::vector<std::ofstream> outClassFSs;
	std::vector<std::ofstream> outClassSimMaxFSs;
	std::vector<std::ofstream> outClassSimAvgFSs;
	std::ofstream outClassMaxFS;
	std::ofstream outnClassSamplesFS;
	std::ofstream outObsSimMaxFS, outObsSimAvgFS;
	uint nClasses;
	bool haveSampleGrid, haveDomainGrid;

	//TESTING GDM SIMILARITY WITH GSTREAMS
	std::vector<igstream> txSrcGSs, txDstGSs;
	
	igstream clsGS, obsGS, inClassSampleGS, inDomainGS;
	std::vector<ogstream> outClassGSs;
	std::vector<ogstream> outClassSimMaxGSs;
	std::vector<ogstream> outClassSimAvgGSs;
	
	ogstream outClassMaxGS;
	ogstream outClassAvgGS;
	ogstream outnClassSamplesGS;
	ogstream outObsSimMaxGS, outObsSimAvgGS;

	ogstream sampleGridGS;

	ogstream outUniquenessGS;

	std::map<uint, uint2> classes;

	//GDM SIMILARITY TO PIXEL
	std::ifstream pixelFS;

	~CBAParams() {

	}
};
#endif
