/*
BDIForRegions.cpp - Calculates total and unique BDI for regions and MBV for classes
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

//TODO replace grid fstreams with gstreams
//TODO tidy up reading and writing vector data
//TODO increase output table floating point precision
//TODO update msgStringStdOut() to include "\n" and align with msgText
//TODO use either msgText() or mggString(). Not both

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <memory>
#include <iomanip>
#include <cmath>
#include "FileUtilsLib.h"
#include "BDIForRegions.h"
#include "BDIUtils.h"

using uint = unsigned int;
using ull = unsigned long long;

#ifndef _WIN32
#	define min std::min
#	define max std::max
#endif

//Calculates BDI for domain or regions and class MBV using class OHA, EHA, grid or probabilities
int CalcBDIForRegions(const char * paramFN) { return CalcBDIForRegions(std::string(paramFN)); }
int CalcBDIForRegions(const std::string & paramFN)
{
	//Filenames and filestreams
	std::string condFN, ehaFN, regFN, clsFN, mbvFN;
	std::string probListFN, ohaTabFN, simTabFN, mbvTabFN, bdiTabFN, BDIDataFN, altOHADataFN;
	std::ifstream condFS, ehaFS, regFS, clsFS;
	std::ofstream bdiTabFS, mbvTabFS, mbvFS;

	//Adding tif probability grid support
	igstream clsGS, tmpGS;

	//Variables
	uint nCols = 0, nRows = 0, nRegions = 1, nClasses = 0;
	size_t nCells = 0;
	double	saRatio = 0.25, mbvAdd = 0.0, mbvFac = 0.0, cellSize = 0.0, areaFactor = 0.0;
	float	noData = -9999.0f;
	uint reg, cls, i, j, k, r, n;
	bool useClassGrid = true, scaleCond = false, multEHACond = false, useBDIData = false, useAltOHAData = false;

	//Get parameters and initialise variables
	msgText("Reading parameters");
	if (!GetParam(paramFN, "BDI", "HAB_FN", condFN)) return msgErrorStr(-1, "HAB_FN", paramFN);
	if (!GetParam(paramFN, "BDI", "EHA_FN", ehaFN)) return msgErrorStr(-1, "EHA_FN", paramFN);
	if (!GetParam(paramFN, "BDI", "VEG_FN", clsFN)) return msgErrorStr(-1, "VEG_FN", paramFN);
	if (!GetParam(paramFN, "BDI", "SAR_Z", saRatio) && !SetParam(paramFN, "BDI", "SAR_Z", saRatio)) return msgErrorStr(-1, "SAR_Z", paramFN);
	if (!GetParam(paramFN, "BDI", "SCALEHAB", scaleCond)) scaleCond = false;
	if (!GetParam(paramFN, "BDI", "MULTEHACOND", multEHACond)) multEHACond = false;
	

	bool useRegions = GetParam(paramFN, "BDI", "REG_FN", regFN);
	bool useOHATab = GetParam(paramFN, "BDI", "OHATABLE", ohaTabFN);
	bool useSimTab = GetParam(paramFN, "BDI", "SIMTABLE", simTabFN) && simTabFN.compare("0") != 0;
	bool useMBVAdd = GetParam(paramFN, "BDI", "MBV_ADD", mbvAdd) && mbvAdd > 0.0;
	bool useMBVFac = GetParam(paramFN, "BDI", "MBV_FAC", mbvFac) && mbvFac > 0.0;
	bool useMBVOHA = !useMBVAdd && !useMBVFac;
	bool createBDITable = GetParam(paramFN, "BDI", "BDITABLE", bdiTabFN);
	bool createMBVTable = GetParam(paramFN, "BDI", "MBVTABLE", mbvTabFN);
	bool createMBVGrid = GetParam(paramFN, "BDI", "MBV_FN", mbvFN);

	if (!GetParam(paramFN, "BDI", "USE_BDI_DATA", useBDIData)) useBDIData = false;
	if (!GetParam(paramFN, "BDI", "USE_ALT_OHA_DATA", useAltOHAData)) useAltOHAData = false;

	if (!GetParam(paramFN, "BDI", "BDI_DATA", BDIDataFN) && useBDIData) return msgErrorStr(-1, "BDI_DATA", paramFN);
	
	if (!GetParam(paramFN, "BDI", "ALT_OHA_DATA", altOHADataFN) && useAltOHAData) return msgErrorStr(-1, "ALT_OHA_DATA", paramFN);


	//Check there is at least one output
	if (!createBDITable && !createMBVTable && !createMBVGrid) 
		return msgErrorStr(-1, "BDITABLE or MBVTABLE or MBV_FN", paramFN);

	//Check condition and EHA grids exist and headers match
	if (!FileExists(GridFN(condFN))) return msgErrorStr(-5, GridFN(condFN));
	if (!FileExists(GridFN(ehaFN))) return msgErrorStr(-5, GridFN(ehaFN));
	if (!ReadGridHeader(condFN, nCols, nRows, cellSize, noData) || !CompareGridHeaders(condFN, ehaFN)) 
		return msgErrorStr(-7, condFN, ehaFN);

	//Check class grid or probability list file exists
	if (clsFN.find(".par") == std::string::npos && FileExists(GridFN(clsFN))) useClassGrid = true;
	else if (GetParam(clsFN, "*", "NPROBGRIDS", nClasses) && nClasses > 0) useClassGrid = false;
	else return msgErrorStr(-99, "Unable to read class grid or probability grid list", clsFN);
	if (useOHATab && (!useClassGrid || useRegions)) return msgErrorStr(-99, "Unable to use OHA table with probability grids or regions");


	//Read regions grid if supplied and count regions and their area including whole region
	nCells = nCols * nRows;
	std::unique_ptr<float[]> regData;
	std::vector<double> regionHa(1, 0.0);
	if (useRegions) {
		float rMax = 0.0f;
		std::cout << "Counting regions\n";
		if (!CompareGridHeaders(condFN, regFN)) return msgErrorStr(-7, condFN, regFN);
		if (!GetGridMax(regFN, rMax) || !(rMax > 0.0f) || rMax == noData) return msgErrorStr(-99, "Unable to get max region from region grid");
		nRegions = int(std::lround(rMax)) + 1;
		regionHa.resize(nRegions, 0.0);
		std::cout << "Number of regions: " << nRegions - 1 << "\n";
		regData = std::make_unique<float[]>(nCells);
		regFS.open(GridFN(regFN), std::ios::in | std::ios::binary);
		if (!regFS.is_open()) return msgErrorStr(-3, regFN);
		regFS.read((char*)regData.get(), nCells * sizeof(float));
		regFS.close();
		for (n = 0; n < nCells; n++) {
			if (regData[n] != noData && regData[n] > 0.0) {
				regionHa[0] += 1.0;
				regionHa[int(std::lround(regData[n]))] += 1.0;
			}
		}
	}

	//Get number of classes from class grid or probability grid filenames from probability grid list
	//Added cond gstream to compare cond and cls grid headers
	tmpGS.open(condFN);
	if (!tmpGS.is_open())
		return msgErrorStr(-3, condFN);

	std::vector<std::string> probFNs;
	if (useClassGrid) {
		msgText("Using class grid");
		float clsMax;
		if (!GetGridMax(clsFN, clsMax) ||
			!CompareGridHeaders(condFN, clsFN)) 
			return msgErrorStr(-99, "Unable to read class grid ", clsFN);
		nClasses = int(std::lround(clsMax));
	}
	else {
		probListFN = clsFN;
		msgText("Using probability grids");
		for (cls = 0; cls < nClasses; cls++) {
			//Added for tif probability grid support
			if (!GetParam(probListFN, "*", "PROBGRID" + toStr(cls + 1), clsFN)) 
				return msgErrorStr(-1, "PROBGRID" + toStr(cls + 1), probListFN);

			clsGS.open(clsFN);
			if (!clsGS.is_open())
				return msgErrorStr(-3, clsFN);
			if (clsGS.compareHeader(tmpGS))
				return msgErrorStr(-7, clsFN, condFN);
			clsGS.close();
			probFNs.push_back(clsFN);
		}
	}
	tmpGS.close();

	//Get class similarities from sim table if using
	std::vector<std::vector<double>> classSim;
	if (useSimTab && !GetClassSimilarities(simTabFN, nClasses, classSim)) return msgErrorStr(-99, "Unable to read similarity table", simTabFN);
	
	//Allocate storage for grid data
	std::unique_ptr<float[]> condData = std::make_unique<float[]>(nCells);
	std::unique_ptr<float[]> ehaData = std::make_unique<float[]>(nCells);
	std::unique_ptr<float[]> clsData = std::make_unique<float[]>(nCells);

	//Read condition and EHA grids
	condFS.open(GridFN(condFN), std::ios::in | std::ios::binary);
	if (!condFS.is_open()) return msgErrorStr(-3, condFN);
	condFS.read((char*)condData.get(), nCells * sizeof(float));
	condFS.close();

	ehaFS.open(GridFN(ehaFN), std::ios::in | std::ios::binary);
	if (!ehaFS.is_open()) return msgErrorStr(-3, ehaFN);
	ehaFS.read((char*)ehaData.get(), nCells * sizeof(float));
	ehaFS.close();

	//Count whole region area if not using region grid
	if (!useRegions) for (i = 0; i < nCells; i++) regionHa[0] += condData[i] == noData ? 0.0 : 1.0;

	//Scale condition to between 0 and 1
	if (scaleCond) {
		float condDen;
		if (!GetGridMaxPow10(condFN, condDen)) return msgErrorStr(-3, condFN);
		if (condDen > 1.0f) {
			msgText(std::string("Scaling condition from between 0 and " + toStr(condDen) + " to between 0 and 1").c_str());
			for (i = 0; i < nCells; i++) condData[i] /= condDen;
		}
	}

	//Multiply EHA by condition 
	//Used by older versions where EHA may exceed condition due to raising cxt to S 
	if (multEHACond) { 
		msgText("Multiplying EHA by condition");
		for (i = 0; i < nCells; i++) ehaData[i] *= condData[i]; 
	}

	//Calculate class condition, EHA and OHA values using class grid
	std::vector<double> classCond(nClasses, 0.0);
	std::vector<std::vector<double>> clsEHAir(nRegions, std::vector<double>(nClasses, 0.0));
	std::vector<std::vector<double>> clsOHAir(nRegions, std::vector<double>(nClasses, 0.0));
	std::vector<std::vector<double>> clsEHAor(nRegions, std::vector<double>(nClasses, 0.0));
	std::vector<std::vector<double>> clsOHAor(nRegions, std::vector<double>(nClasses, 0.0));
	
	if (useBDIData) {
		msgString("Reading class condition, EHA and OHA from " + BDIDataFN + "_*.bin\n");
		if (!readVectorFromFile<double>(BDIDataFN + "_clscond.bin", classCond)) return msgErrorStr(-3, BDIDataFN + "_clscond.bin");
		if (!read2dVectorFromFile<double>(BDIDataFN + "_ehair.bin", clsEHAir)) return msgErrorStr(-3, BDIDataFN + "_ehair.bin");
		if (!read2dVectorFromFile<double>(BDIDataFN + "_ohair.bin", clsOHAir)) return msgErrorStr(-3, BDIDataFN + "_ohair.bin");
		if (!read2dVectorFromFile<double>(BDIDataFN + "_ehaor.bin", clsEHAor)) return msgErrorStr(-3, BDIDataFN + "_ehaor.bin");
		if (!read2dVectorFromFile<double>(BDIDataFN + "_ohaor.bin", clsOHAor)) return msgErrorStr(-3, BDIDataFN + "_ohaor.bin");
	}
	else {
		double cellEHA = 0.0;
		if (useClassGrid) {
			msgText("Calculating class condition, EHA and OHA using class grid");
			clsFS.open(GridFN(clsFN), std::ios::in | std::ios::binary);
			if (!clsFS.is_open()) return msgErrorStr(-3, clsFN);
			clsFS.read((char*)clsData.get(), nCells * sizeof(float));
			clsFS.close();
			for (i = 0; i < nCells; i++) {
				if (clsData[i] != noData && condData[i] != noData && ehaData[i] != noData && int(std::lround(clsData[i])) > 0) {
					
					//Moved regData checks up to ensure EHA/OHA[0] sums only include pixels in regions 
					if (useRegions && (regData[i] == noData || int(std::lround(regData[i])) < 1)) continue;
					
					cls = uint(std::lround(clsData[i])) - 1;
					cellEHA = double(ehaData[i]);
					classCond[cls] += double(condData[i]);
					clsEHAir[0][cls] += cellEHA;
					clsOHAir[0][cls] += 1.0;

					//Calculate class EHA and OHA values for region at cell i
					if (useRegions) { // && regData[i] != noData && int(std::lround(regData[i])) > 0) {//No data and > 0 now checked above
						reg = uint(std::lround(regData[i]));
						clsEHAir[reg][cls] += cellEHA;
						clsOHAir[reg][cls] += 1.0;
						for (r = 1; r < nRegions; r++) {
							if (r != reg) {
								clsEHAor[r][cls] += cellEHA;
								clsOHAor[r][cls] += 1.0;
							}
						}
					}
				}
			}
		}

		//Calculate class condition, EHA and OHA values using class probability grids
		else {
			msgText("Calculating class condition, EHA and OHA using probability grids");
			double cellProb = 0.0;
			for (cls = 0; cls < nClasses; cls++) {
				msgProgress("Percent complete: ", int(100 * cls / nClasses));
				//Added tif probability grid support
				clsGS.open(probFNs[cls]);
				if (!clsGS.is_open()) return msgErrorStr(-3, probFNs[cls]);
				clsGS.read((char*)clsData.get(), nCells * sizeof(float));
				clsGS.close();

				for (i = 0; i < nCells; i++) {
					if (clsData[i] > 0.0f && clsData[i] != noData && condData[i] != noData && ehaData[i] != noData) {
						
						//Moved regData checks up to ensure EHA/OHA[0] sums only include pixels in regions
						// TODO this breaks benefit grids where no original extent is within regions
						//Can't use regions with incomplete coverage and derive full MBV coverage
						if (useRegions && (regData[i] == noData || int(std::lround(regData[i])) < 1)) continue;
					
						cellProb = double(clsData[i]);
						cellEHA = cellProb * double(ehaData[i]);
						classCond[cls] += cellProb * double(condData[i]);
						clsEHAir[0][cls] += cellEHA;
						clsOHAir[0][cls] += cellProb;

						//Calculate class EHA and OHA values for region at cell i
						if (useRegions) { // && regData[i] != noData && int(std::lround(regData[i])) > 0) {//No data and > 0 now checked above
							reg = uint(std::lround(regData[i]));
							clsEHAir[reg][cls] += cellEHA;
							clsOHAir[reg][cls] += cellProb;
							for (r = 1; r < nRegions; r++) {
								if (r != reg) {
									clsEHAor[r][cls] += cellEHA;
									clsOHAor[r][cls] += cellProb;
								}
							}
						}
					}
				}
			}
			msgText("\rPercent complete: 100");
		}
		
		msgString("Writing class condition, EHA and OHA to " + BDIDataFN + "_*.bin\n");
		if (!writeVectorToFile(BDIDataFN + "_clscond.bin", classCond)) return msgErrorStr(-4, BDIDataFN + "_clscond.bin");
		if (!write2dVectorToFile(BDIDataFN + "_ehair.bin", clsEHAir)) return msgErrorStr(-4, BDIDataFN + "_ehair.bin");
		if (!write2dVectorToFile(BDIDataFN + "_ohair.bin", clsOHAir)) return msgErrorStr(-4, BDIDataFN + "_ohair.bin");
		if (!write2dVectorToFile(BDIDataFN + "_ehaor.bin", clsEHAor)) return msgErrorStr(-4, BDIDataFN + "_ehaor.bin");
		if (!write2dVectorToFile(BDIDataFN + "_ohaor.bin", clsOHAor)) return msgErrorStr(-4, BDIDataFN + "_ohaor.bin");
	}

	//Sum region condition and EHA values for all regions
	std::vector<double> regionCond(nRegions, 0.0); //Condition for each region;
	std::vector<double> regionEHA(nRegions, 0.0); //EHA for each region
	
	for (i = 0; i < nCells; i++) {
		if (condData[i] != noData && ehaData[i] != noData) {
			regionCond[0] += double(condData[i]);
			regionEHA[0] += double(ehaData[i]);
			if (useRegions && regData[i] != noData && int(std::lround(regData[i])) > 0) {
				reg = uint(std::lround(regData[i]));
				regionCond[reg] += double(condData[i]);
				regionEHA[reg] += double(ehaData[i]);
			}
		}
	}

	//Read alternate OHA data if using
	if (useAltOHAData) {
		msgString("Reading alternate OHA data from " + altOHADataFN + "_*.bin\n");
		if (!read2dVectorFromFile<double>(altOHADataFN + "_ohair.bin", clsOHAir)) return msgErrorStr(-3, altOHADataFN + "_ohair.bin");
		if (!read2dVectorFromFile<double>(altOHADataFN + "_ohaor.bin", clsOHAor)) return msgErrorStr(-3, altOHADataFN + "_ohaor.bin");
	}

	//Reset condition, EHA and region grid data pointers to remove from memory
	condData.reset();
	ehaData.reset();
	regData.reset();
	
	//Update class OHA from OHA table if using
	if (useOHATab && !GetClassOHAFromTable(ohaTabFN, nClasses, clsOHAir)) return msgErrorStr(-99, "Unable to read class OHA values from", ohaTabFN);

	//Apply areaFactor to convert values to hectares and calculate increased EHA for class MBV
	areaFactor = pow(cellSize, 2) * (cellSize > 1.0f ? 0.0001L : 1000000.0L);
	std::vector<double> clsMBVEHA(nClasses, 0.0);

	//Divide region conditon and EHA sums by count to get averages before scaling count by area factor
	for (r = 0; r < nRegions; r++) regionCond[r] /= regionHa[r];
	for (r = 0; r < nRegions; r++) regionEHA[r] /= regionHa[r];
	for (r = 0; r < nRegions; r++) regionHa[r] *= areaFactor;

	for (n = 0; n < nClasses; n++) {
		classCond[n] *= areaFactor;
		for (r = 0; r < nRegions; r++) {
			clsEHAir[r][n] *= areaFactor;
			clsOHAir[r][n] *= useOHATab ? 1.0 : areaFactor;
			clsEHAor[r][n] *= areaFactor;
			clsOHAor[r][n] *= areaFactor;
		}
		if (useMBVAdd) clsMBVEHA[n] = clsEHAir[0][n] + mbvAdd;
		//else if (useMBVFac) clsMBVEHA[n] = clsEHAir[0][n] * mbvFac; Add % of OHA rather than factor of EHA
		else if (useMBVFac) clsMBVEHA[n] = clsEHAir[0][n] + (clsOHAir[0][n] * mbvFac);
		else clsMBVEHA[n] = clsOHAir[0][n];
	}

	//TODO this could be split into separate functions. 

	//Storage for BDI and MBV values
	//ir = in region, nr = not in region, pr = pristine in region, pn = pristine not in region
	msgText((std::string("Calculating BDI and MBV ") +
		(useSimTab ? "using similarities and " : "") +
		(useMBVAdd ? "MBV addition" : useMBVFac ? "MBV factor" : "OHA")).c_str());
	double irBDINum = 0.0;
	double nrBDINum = 0.0;
	double prBDINum = 0.0;
	double pnBDINum = 0.0;
	double BDIDen = 0.0;
	std::vector<double> irBDI(nRegions, 0.0);
	std::vector<double> nrBDI(nRegions, 0.0);
	std::vector<double> prBDI(nRegions, 0.0);
	std::vector<double> pnBDI(nRegions, 0.0);

	std::vector<double> classBDI(nClasses, 0.0);
	std::vector<double> classMBV(nClasses, 0.0);

	//Calculate region BDI and class MBV values without similarities
	if (!useSimTab) {
		//TODO Need to check these calcs
		//TODO Calc BDI for each class before summing for region
		return msgErrorStr(-99, "Calculating region BDI and class MBV values without similarities not currently implemented");
		/*
		for (r = 0; r < nRegions; r++) {
			irBDINum = 0.0;
			nrBDINum = 0.0;
			prBDINum = 0.0;
			pnBDINum = 0.0;
			BDIDen = 0.0;
			for (i = 0; i < nClasses; i++) {
				irBDINum += clsEHAir[r][i];
				nrBDINum += clsEHAor[r][i];
				prBDINum += clsOHAir[r][i];
				pnBDINum += clsOHAor[r][i];
				BDIDen += clsOHAir[0][i];
			}
			if (BDIDen > 0.0) {
				irBDI[r] = pow(irBDINum / BDIDen, saRatio);
				nrBDI[r] = pow(nrBDINum / BDIDen, saRatio);
				prBDI[r] = pow(prBDINum / BDIDen, saRatio);
				pnBDI[r] = pow(pnBDINum / BDIDen, saRatio);
			}
		}
			
		//Calculte MBV for each class in region 0
		for (n = 0; n < nClasses; n++) {
			irBDINum = 0.0;
			BDIDen = 0.0;
			for (i = 0; i < nClasses; i++) {
				irBDINum += (i == n) ? clsMBVEHA[i] : clsEHAir[0][i];
				BDIDen += clsOHAir[0][i];
			}
			if (BDIDen > 0.0) classMBV[n] = pow(irBDINum / BDIDen, saRatio) - irBDI[0];
		}
		*/
		
	}

	//Calculate region BDI and class MBV values with similarities
	else {
		double sumEHAjir = 0.0;
		double sumEHAjor = 0.0;
		double sumOHAjir = 0.0;
		double sumOHAjor = 0.0;
		double sumOHAjar = 0.0;

		double clsBDINum = 0.0, clsBDIDen = 0.0;

		for (r = 0; r < nRegions; r++) {
			irBDINum = 0.0;
			nrBDINum = 0.0;
			prBDINum = 0.0;
			pnBDINum = 0.0;
			BDIDen = 0.0;

			//Calculate BDI for region r including not region (Nr), pristine region (Pr) and pristine not region (Pn)
			for (i = 0; i < nClasses; i++) {
				if (clsOHAir[0][i] <= 0.0) continue;

				//Zero for each i
				clsBDINum = 0.0;
				clsBDIDen = 0.0;

				for (k = 0; k < 11; k++) {
					sumEHAjir = 0.0;
					sumEHAjor = 0.0;
					sumOHAjir = 0.0;
					sumOHAjor = 0.0;
					sumOHAjar = 0.0;
					for (j = 0; j < nClasses; j++) {
						if (classSim[i][j] >= (double(k) * 0.1)) {
							sumEHAjir += clsEHAir[r][j];
							sumEHAjor += clsEHAor[r][j];
							sumOHAjir += clsOHAir[r][j];
							sumOHAjor += clsOHAor[r][j];
							sumOHAjar += clsOHAir[0][j];
						}
					}// for j
					if (sumOHAjar > 0.0) {
						sumEHAjir = min(sumEHAjir, sumOHAjar);//Clamp EHA to OHA, Needed?
						clsBDINum += sumEHAjir;
						clsBDIDen += sumOHAjar;
						irBDINum += clsOHAir[0][i] * pow(sumEHAjir / sumOHAjar, saRatio) / sumOHAjar;
						nrBDINum += clsOHAir[0][i] * pow(sumEHAjor / sumOHAjar, saRatio) / sumOHAjar;
						prBDINum += clsOHAir[0][i] * pow(sumOHAjir / sumOHAjar, saRatio) / sumOHAjar;
						pnBDINum += clsOHAir[0][i] * pow(sumOHAjor / sumOHAjar, saRatio) / sumOHAjar;
						BDIDen += clsOHAir[0][i] / sumOHAjar;
					}
				}// for k
				if (r == 0) {
					classBDI[i] = pow(clsBDINum / clsBDIDen, saRatio);
				}
			}// for i
			if (BDIDen > 0.0) {
				irBDI[r] = irBDINum / BDIDen;
				nrBDI[r] = nrBDINum / BDIDen;
				prBDI[r] = prBDINum / BDIDen;
				pnBDI[r] = pnBDINum / BDIDen;
			}
		}//for r

		//Calculte MBV for each class cls in region 0 using similarities
		for (n = 0; n < nClasses; n++) {
			irBDINum = 0.0;
			BDIDen = 0.0;
			for (i = 0; i < nClasses; i++) {
				if (clsOHAir[0][i] <= 0.0) continue;
				for (k = 0; k < 11; k++) {
					sumEHAjir = 0.0;
					sumOHAjar = 0.0;
					for (j = 0; j < nClasses; j++) {
						if (classSim[i][j] >= (double(k) * 0.1)) {
							sumEHAjir += (j == n) ? clsMBVEHA[j] : clsEHAir[0][j];
							sumOHAjar += clsOHAir[0][j];
						}
					}// for j
					if (sumOHAjar > 0.0) {
						sumEHAjir = min(sumEHAjir, sumOHAjar);//Clamp EHA to OHA, Needed?
						irBDINum += clsOHAir[0][i] * pow(sumEHAjir / sumOHAjar, saRatio) / sumOHAjar;
						BDIDen += clsOHAir[0][i] / sumOHAjar;
					}
				}// for k
			}// for i
			if (BDIDen > 0.0) classMBV[n] = (irBDINum / BDIDen) - irBDI[0];
		}// for n
	}

	//Calculate total and unique current and pristine diversity for each region
	std::vector<double> totalBDI(nRegions, 0.0);
	std::vector<double> totalPristBDI(nRegions, 0.0);
	std::vector<double> uniqueBDI(nRegions, 0.0);
	std::vector<double> uniquePristBDI(nRegions, 0.0);
	for (r = 0; r < nRegions; r++) {
		totalBDI[r] = irBDI[r] * 100.0;
		totalPristBDI[r] = prBDI[r] * 100.0;
		uniqueBDI[r] = (irBDI[0] - nrBDI[r]) * 100.0;
		uniquePristBDI[r] = (prBDI[0] - pnBDI[r]) * 100.0;
	}

	//Write BDI table
	if (createBDITable) {
		msgText("Creating BDI table");
		bdiTabFS.open(bdiTabFN, std::ios::out);
		if (!bdiTabFS.is_open()) return msgErrorStr(-4, bdiTabFN);
		bdiTabFS << "Region, Area (ha), Condition (Avg), EHA (Avg), Total diversity (Pre),Total diversity (Cur),Unique diversity (Pre),Unique diversity (Cur)\n";
		for (r = 0; r < nRegions; r++) {
			SetParam(paramFN, "BDI", "BDI" + toStr(r), irBDI[r]);
			bdiTabFS << r << "," << regionHa[r] << "," << regionCond[r] << "," << regionEHA[r] << "," << totalPristBDI[r] << "," << totalBDI[r] << "," << uniquePristBDI[r] << "," << uniqueBDI[r] << "\n";
		}
		bdiTabFS.close();
	}

	//Write MBV table
	if (createMBVTable) {
		msgText("Creating MBV table");
		mbvTabFS.open(mbvTabFN, std::ios::out);
		if (!mbvTabFS.is_open()) return msgErrorStr(-4, mbvTabFN);
		mbvTabFS << "CLASS,OHA (Ha),EHA (Ha),CONDITION (Ha),BDI,MBV\n";
		for (cls = 0; cls < nClasses; cls++)
			mbvTabFS << cls + 1 << "," << clsOHAir[0][cls] << "," << clsEHAir[0][cls] << "," << classCond[cls] << "," << classBDI[cls] << "," << classMBV[cls] << "\n";
		mbvTabFS.close();
	}

	//Create MBV grid
	if (createMBVGrid) {
		std::unique_ptr<float[]> mbvData = std::make_unique<float[]>(nCells);
		mbvFS.open(GridFN(mbvFN), std::ios::out | std::ios::binary);
		if (!mbvFS.is_open()) return msgErrorStr(-4, GridFN(mbvFN));
		
		//Create MBV grid using class grid
		if (useClassGrid) {
			msgText("Creating MBV grid from class grid");
			for (i = 0; i < nCells; i++)
				mbvData[i] = clsData[i] == noData ?	noData : float(classMBV[uint(std::lround(clsData[i])) - 1]);
		}
		
		//Create MBV grid using class probability grids
		else {
			msgText("Creating MBV grid from probability grids");
			for (i = 0; i < nCells; i++) mbvData[i] = 0.0f;
			for (cls = 0; cls < nClasses; cls++) {
				msgProgress("Percent complete: ", int(100 * cls / nClasses));
				//Added tif probability grid support
				clsGS.open(probFNs[cls]);
				clsGS.read((char*)clsData.get(), nCells * sizeof(float));
				clsGS.close();
				for (i = 0; i < nCells; i++) {
					if (clsData[i] == noData) mbvData[i] = noData;
					else mbvData[i] += float(double(clsData[i]) * (classMBV[cls]));
				}
			}
			msgText("\rPercent complete: 100");
		}
		
		//Write MBV grid to file
		mbvFS.write((const char*)mbvData.get(), nCells * sizeof(float));
		mbvFS.close();
		CopyFS(HeaderFN(condFN), HeaderFN(mbvFN), true);
	}

	//Complete!
	msgText("CalcBDIForRegions() Complete!");
	return 0;
}

