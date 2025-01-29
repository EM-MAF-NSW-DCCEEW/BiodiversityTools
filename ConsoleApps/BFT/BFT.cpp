/*
BFT.cpp - Biodiversity Foreacsting Tool 
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

/* TODO
TEST - Use provided Prm Min/Max
Add region name lookup table
*/


#include <fstream>
#include <iostream>
#include <cmath>
#include "FileUtilsLib.h"
#include "SpatialContextLib.h"
#include "BDIForRegions.h"


int RunBFT(std::string paramFN);

int main(int argc, char *argv[])
{
	//Options
	//Pause after use
	if (argc < 2) {
		std::cout << "\nUsage: " << argv[0] << " ParameterFile.par \n"
			<< "\nIf ParameterFile.par exists it is processed, otherwise it is created as a template.\n\n";
		return -1;
	}
	std::fstream testFS;
	std::string argv1(argv[1]);
	testFS.open(argv1, std::ios::in);
	if (!testFS.is_open()) {
		std::cout << "\nThe parameter file " << argv[1] << " does not exist and is being created as a template\n\n";
		if (
			//Read parameters
			!SetParam(argv1, "INPUTS", "HAB_FN", "Filename. Habitat input grid file name without extension.") ||
			!SetParam(argv1, "INPUTS", "HABMAX", "Floating point number. Maximum possible habitat value. Optional, default = smallest power of 10 greater than maximum habitat grid value.") || 
			!SetParam(argv1, "INPUTS", "PRM_FN", "Filename. Permeability input grid file name without extension. Optional, default calculated from habitat grid using DISTMIN and DISTMAX") ||
			!SetParam(argv1, "INPUTS", "EHA_FN", "Filename. EHA input grid filename without extension. Optional, default calculated") ||
			!SetParam(argv1, "INPUTS", "MBV_FN", "Filename. MBV input grid filename without extension. Optional, default calculated") ||

			!SetParam(argv1, "WINDOW", "SRCSIZE", "Integer. Source window size in cells. Optional, default = use existing lookup and translate tables.") ||
			!SetParam(argv1, "WINDOW", "DSTSIZE", "Integer. Destination window size in petals. Optional, default = use existing lookup and translate tables.") ||
			!SetParam(argv1, "WINDOW", "ZONERATIO", "Floating point number. The area ratio increase of each petal outwards. Must be >= 1.0. Optional, default = use existing lookup and translate tables.") ||
			!SetParam(argv1, "WINDOW", "SAVETMPFILES", "true/false. Determines whether intermediate files are created. Optional, default = false.") ||
			!SetParam(argv1, "WINDOW", "LOOKUP_FN", "Filename. Petals look up file name with extension. Source and destination window sizes are derived from this. Optional, default = use window parameters.") ||
			!SetParam(argv1, "WINDOW", "TRANSLATE_FN", "Filename. Petals translate file name with extension. Source and destination window sizes are derived from this. Optional, default = use window parameters.") ||

			!SetParam(argv1, "ANALYSIS", "WORKDIR", "Path to an existing folder where outputs are written to. Optional, default = parameter file's folder if resolvable.") ||
			!SetParam(argv1, "ANALYSIS", "DISTMIN", "Floating point number. Minimum (degraded) average dispersal distance. Optional, either DISTMIN/MAX or PERMMIN/MAX required.") ||
			!SetParam(argv1, "ANALYSIS", "DISTMAX", "Floating point number. Maximum (intact) average dispersal distance. Optional, either DISTMIN/MAX or PERMMIN/MAX required.") ||
			!SetParam(argv1, "ANALYSIS", "PERMMIN", "Floating point number. Minimum (degraded) average dispersal distance. Optional, either DISTMIN/MAX or PERMMIN/MAX required.") ||
			!SetParam(argv1, "ANALYSIS", "PERMMAX", "Floating point number. Maximum (intact) average dispersal distance. Optional, either DISTMIN/MAX or PERMMIN/MAX required.") ||
			!SetParam(argv1, "ANALYSIS", "MULTFOCAL", "true/false, Whether to multiply connectivity by the focal cell habitat value. Optional, default = false.") ||
			!SetParam(argv1, "ANALYSIS", "FOCALPOWER", "Floating point number. Exponent applied to the focal cell habitat value. Optional, default = 1.0.") ||
			!SetParam(argv1, "ANALYSIS", "SUMPOWER", "Floating point number. Exponent applied to EHA. Optional, default = 1.0.") ||

			//BDI/MBV Assessment parameters
			!SetParam(argv1, "ASSESSMENT", "REG_FN", "Filename. Region input grid file name without extension. Grid values from 1 to n for n regions. Optional. default = not used.") ||
			!SetParam(argv1, "ASSESSMENT", "VEG_FN", "Filename. Veg class input grid file name without extension or .par file containing NPROBGRIDS and PROBGRIDn for n = 1 to NPROBGRIDS.") ||
			!SetParam(argv1, "ASSESSMENT", "SIMTABLE", "Filename. Class similarity input table file name with .csv extension. Optional, default = not used.") ||
			!SetParam(argv1, "ASSESSMENT", "OHATABLE", "Filename. OHA table with .csv extension. Optional, can only be used with class grid and no regions.") ||
			!SetParam(argv1, "ASSESSMENT", "SCALEHAB", "true/false. Whether to scale the habitat grid to between 0 and 1. Optional, default = false.") ||
			!SetParam(argv1, "ASSESSMENT", "MULTEHACOND", "true/false. Whether to multiply EHA by condition. Optional, default = false.") ||
			!SetParam(argv1, "ASSESSMENT", "SAR_Z", "Floating point number. Species area curve exponent. Optional, default = 2.5.") ||
			!SetParam(argv1, "ASSESSMENT", "MBV_ADD", "Floating point number. Area (ha) to add to EHA for each class when calculating MBV. Optional, default = MBV_FAC or use OHA.") ||
			!SetParam(argv1, "ASSESSMENT", "MBV_FAC", "Floating point number. Percent of OHA area less than 1.0 to apply to EHA for each class when calculating MBV. Optional, default = MBV_ADD or use OHA.") ||
			!SetParam(argv1, "ASSESSMENT", "BDITABLE", "Filename. BDI table output file name with .csv extension. Optional, default = not used.") ||
			!SetParam(argv1, "ASSESSMENT", "MBVTABLE", "Filename. MBV table output file name with .csv extension. Optional, default = not used.") ||
			!SetParam(argv1, "ASSESSMENT", "MBV_FN", "Filename. MBV output grid file name without extension. Optional, default = not used.") ||

			//!SetParam(argv1, "ASSESSMENT", "USE_BDI_DATA", "TODO") ||
			!SetParam(argv1, "ASSESSMENT", "BDI_DATA", "Filename. Use existing BDI data files with specified filename prefix for BDI and MBV calcs if avaiable. Optional, default = not used.") ||

			//!SetParam(argv1, "ASSESSMENT", "USE_ALT_OHA_DATA", "TODO") ||
			!SetParam(argv1, "ASSESSMENT", "ALT_OHA_DATA", "Filename. Use alternate OHA data files for BDI and MBV calcs if avaiable. Optional, default = not used.") ||


			!SetParam(argv1, "ASSESSMENT", "BENEFITS", "Comma seperated text values. List of benefit grids to generate each with a matching named parameter section as per following. Optional, default = not used.") ||
			!SetParam(argv1, "ASSESSMENT", "#NOTE", "If 'BENEFITS=conserve, improve, restore' then [conserve], [improve] and [restore] parameter sections containing the following [BENEFIT] parameters are required") ||

			!SetParam(argv1, "BENEFIT", "BEN_FN", "Filename. Benefit output grid filename without extension. Optional, default = benefit parameter section name") ||
			!SetParam(argv1, "BENEFIT", "MBV_FN", "Filename. Alternate MBV grid filename without extension. Optional, default = MBV_FN from [ASSESSMENT]") ||
			!SetParam(argv1, "BENEFIT", "ALTGRIDS", "true/false. ALTHAB and ALTPRM are alternate input grid filenames or constant values if false.") ||
			!SetParam(argv1, "BENEFIT", "ALTHAB", "Filename or floating point number. Alternate habitat input grid file name without extension when ALTGRIDS=true or constant value when ALTGRIDS=false.") ||
			!SetParam(argv1, "BENEFIT", "ALTPRM", "Filename or floating point number. Alternate permeability input grid file name without extensionwhen ALTGRIDS=true or constant value when ALTGRIDS=false.") ||
			!SetParam(argv1, "BENEFIT", "MULTFOCAL", "true/false, Whether to multiply connectivity by the focal cell habitat value. Optional, default = MULTFOCAL from [ANALYSIS]") ||
			!SetParam(argv1, "BENEFIT", "FOCALPOWER", "Floating point number. Exponent applied to the focal cell habitat value. Optional, default = SUMPOWER from [ANALYSIS]") ||
			!SetParam(argv1, "BENEFIT", "SUMPOWER", "Floating point number. Exponent applied to EHA. Optional, default = MBV_FN from [ANALYSIS]")
			)
			return msgErrorStr(-99, "Unable to create parameter file template:", argv1);
		return -2;
	}
	testFS.close();
	int returnValue = RunBFT(argv1);
	if (returnValue != 0) std::cout << "Error code: " << returnValue << "\n";
	else std::cout << "Complete!" << "\n";
	return returnValue;
}


int RunBFT(std::string paramFN) {

	std::string paramDir("");
	std::string workDir("");
	
	bool haveHabMax = false;
	bool createPrmGrid = false;
	bool createEHAGrid = false;
	bool doBDIAndMBV = false;
	bool multFocal = false;
	bool scaleCond = false;
	bool multEHACond = false;
	bool doBenefits = false;

	float habMax = 0.0f;
	float distMin = 0.0f, distMax = 0.0f;
	double prmMin = 0.0, prmMax = 0.0;
	float focalPower, sumPower;
	
	int result;

	std::string habFN, prmFN, ehaFN;
	std::string regFN, clsFN, mbvFN, probFN;
	std::string probListFN, ohaTabFN, simTabFN, mbvTabFN, bdiTabFN, tmpDir;
	std::string BDIDataFN = "bdi";
	std::string altOHADataFN;
	//std::ifstream condFS, ehaFS, regFS;
	//std::ofstream bdiTabFS, mbvTabFS, mbvFS;

	double	saRatio = 0.25, mbvAdd = 0.0, mbvFac = 0.0, cellSize = 0.0, areaFactor = 0.0;

	//Read parameters
	if (!GetParam(paramFN, "ANALYSIS", "WORKDIR", workDir) && (!FileExists(paramFN) || !GetFilesFolder(paramFN, workDir)))
		return msgErrorStr(-1, "file's directory or WORKDIR", paramFN);

	if (!GetParam(paramFN, "INPUTS", "HAB_FN", habFN)) return msgErrorStr(-1, "HAB_FN", paramFN);
	if (GetParam(paramFN, "INPUTS", "HABMAX", habMax)) haveHabMax = true;

	if (!(GetParam(paramFN, "INPUTS", "PRM_FN", prmFN) && FileExists(GridFN(prmFN)))) createPrmGrid = true;
	if (!(GetParam(paramFN, "INPUTS", "EHA_FN", ehaFN) && FileExists(GridFN(ehaFN)))) createEHAGrid = true;
	if (!(GetParam(paramFN, "INPUTS", "MBV_FN", mbvFN) && FileExists(GridFN(mbvFN)))) doBDIAndMBV = true;

	if (createPrmGrid && 
		!((GetParam(paramFN, "ANALYSIS", "DISTMIN", distMin) && GetParam(paramFN, "ANALYSIS", "DISTMAX", distMax)) ||
		((GetParam(paramFN, "ANALYSIS", "PERMMIN", prmMin) && GetParam(paramFN, "ANALYSIS", "PERMMAX", prmMax))))) 
		return msgErrorStr(-1, "DISTMIN/MAX or PERMMIN/MAX", paramFN);

	if (!GetParam(paramFN, "ANALYSIS", "MULTFOCAL", multFocal)) multFocal = false;
	if (!GetParam(paramFN, "ANALYSIS", "FOCALPOWER", focalPower)) focalPower = 1.0f;
	if (!GetParam(paramFN, "ANALYSIS", "SUMPOWER", sumPower)) sumPower = 1.0f;

	//BDI/MBV Assessment parameters
	if (!GetParam(paramFN, "ASSESSMENT", "VEG_FN", clsFN)) return msgErrorStr(-1, "VEG_FN", paramFN);
	if (!GetParam(paramFN, "ASSESSMENT", "SAR_Z", saRatio) && !SetParam(paramFN, "ASSESSMENT", "SAR_Z", saRatio)) return msgErrorStr(-1, "SAR_Z", paramFN);
	if (!GetParam(paramFN, "ASSESSMENT", "SCALEHAB", scaleCond)) scaleCond = false;
	if (!GetParam(paramFN, "ASSESSMENT", "MULTEHACOND", multEHACond)) multEHACond = false;

	bool useRegions = GetParam(paramFN, "ASSESSMENT", "REG_FN", regFN);
	bool useOHATab = GetParam(paramFN, "ASSESSMENT", "OHATABLE", ohaTabFN);
	bool useSimTab = GetParam(paramFN, "ASSESSMENT", "SIMTABLE", simTabFN) && simTabFN.compare("0") != 0;
	bool useMBVAdd = GetParam(paramFN, "ASSESSMENT", "MBV_ADD", mbvAdd) && mbvAdd > 0.0;
	bool useMBVFac = GetParam(paramFN, "ASSESSMENT", "MBV_FAC", mbvFac) && mbvFac > 0.0;
	
	bool createBDITable = doBDIAndMBV && GetParam(paramFN, "ASSESSMENT", "BDITABLE", bdiTabFN);
	bool createMBVTable = doBDIAndMBV && GetParam(paramFN, "ASSESSMENT", "MBVTABLE", mbvTabFN);
	bool createMBVGrid = doBDIAndMBV && GetParam(paramFN, "ASSESSMENT", "MBV_FN", mbvFN);

	bool useBDIData = doBDIAndMBV && GetParam(paramFN, "ASSESSMENT", "BDI_DATA", BDIDataFN);
	bool useAltOHAData = doBDIAndMBV && GetParam(paramFN, "ASSESSMENT", "ALT_OHA_DATA", altOHADataFN);

	if (createBDITable && !GetFilesFolder(bdiTabFN, tmpDir)) bdiTabFN = workDir + pathSep + bdiTabFN;
	if (createMBVTable && !GetFilesFolder(mbvTabFN, tmpDir)) mbvTabFN = workDir + pathSep + mbvTabFN;
	if (createMBVGrid && !GetFilesFolder(mbvFN, tmpDir)) mbvFN = workDir + pathSep + mbvFN;
	
	BDIDataFN = workDir + pathSep + BDIDataFN;

	std::vector<std::string> benefits;
	if (GetParam(paramFN, "ASSESSMENT", "BENEFITS", benefits)) doBenefits = true;

	std::cout << "Creating work directory\n";
	if (!(fileSys::is_directory(fileSys::path(workDir))) && !(fileSys::create_directories(fileSys::path(workDir))))
		return msgErrorStr(-99, "Unable to find or create folder:", workDir);

	std::string lookupFN(""), translateFN("");
	if (createEHAGrid || doBenefits) {
		//Get petal files
		if (!(ParamSectionExists(paramFN, "WINDOW") && GetPetalFiles(paramFN, "WINDOW", workDir, lookupFN, translateFN) == 0))
			return msgErrorStr(-99, "Unable to get petals from", paramFN);
		
		//Get maximum habitat value and cell size
		if (!haveHabMax && !GetGridMaxPow10(habFN, habMax)) return msgErrorStr(-99, "Unable to get habitat max from", habFN);
		if (!GetHeader(habFN, "cellsize", cellSize)) return msgErrorStr(-5, habFN);

		//Create permeability grid
		if (createPrmGrid) {
			prmFN = workDir + pathSep + "prm";
			if (distMin > 0.0f) prmMin = std::exp(-(cellSize < 1.0 ? cellSize * 100000 : cellSize) / distMin);
			if (distMax > 0.0f) prmMax = std::exp(-(cellSize < 1.0 ? cellSize * 100000 : cellSize) / distMax);
			std::cout << "Creating permeability grid " << prmFN << " with values from " << prmMin << " to " << prmMax << "\n";
			if (!RescaleGrid(habFN, prmFN, 0.0f, habMax, float(prmMin), float(prmMax))) return msgErrorStr(-13, prmFN);
			if (!SetParam(paramFN, "INPUTS", "PRM_FN", prmFN)) return msgErrorStr(-14, "PRM_FN", paramFN);
		}
	}

	//Run EHA
	if (createEHAGrid) {
		ehaFN = std::string(workDir + pathSep + "eha");
		std::string ehaParamFN(ehaFN + ".par");

		std::cout << "Creating EHA parameter file: " << ehaParamFN << "\n";
		if (!SetParam(ehaParamFN, "INPUTS", "HAB_FN", habFN) ||
			!SetParam(ehaParamFN, "INPUTS", "PRM_FN", prmFN) ||
			!SetParam(ehaParamFN, "OUTPUTS", "EHA_FN", ehaFN) ||
			!SetParam(ehaParamFN, "WINDOW", "LOOKUP_FN", lookupFN) ||
			!SetParam(ehaParamFN, "WINDOW", "TRANSLATE_FN", translateFN) ||
			!SetParam(ehaParamFN, "ANALYSIS", "MULTFOCAL", multFocal) ||
			!SetParam(ehaParamFN, "ANALYSIS", "FOCALPOWER", focalPower) ||
			!SetParam(ehaParamFN, "ANALYSIS", "SUMPOWER", sumPower))
			return msgErrorStr(-14, "parameters", ehaParamFN);
		
		std::cout << "Processing EHA parameter file: " << ehaParamFN << "\n";
		result = EHA_CBA(ehaParamFN.c_str());
		if (result != 0) return msgErrorStr(-15, "EHA_CBA()", toStr(result));
		if (!SetParam(paramFN, "ASSESSMENT", "EHA_FN", ehaFN)) return msgErrorStr(-14, "EHA_FN", paramFN);
	}

	//Run BDI and MBV
	if (doBDIAndMBV) {
		std::string bdiParamFN(workDir + pathSep + "bdi.par");

		std::cout << "Creating BDI parameter file: " << bdiParamFN << "\n";
		if (!SetParam(bdiParamFN, "BDI", "HAB_FN", habFN) ||
			!SetParam(bdiParamFN, "BDI", "EHA_FN", ehaFN) ||
			!SetParam(bdiParamFN, "BDI", "VEG_FN", clsFN) ||
			(useRegions && !SetParam(bdiParamFN, "BDI", "REG_FN", regFN)) ||
			(useOHATab && !SetParam(bdiParamFN, "BDI", "OHATABLE", ohaTabFN)) ||
			(useSimTab && !SetParam(bdiParamFN, "BDI", "SIMTABLE", simTabFN)) ||
			
			(createBDITable && !SetParam(bdiParamFN, "BDI", "BDITABLE", bdiTabFN)) ||
			(createMBVTable && !SetParam(bdiParamFN, "BDI", "MBVTABLE", mbvTabFN)) ||
			(createMBVGrid && !SetParam(bdiParamFN, "BDI", "MBV_FN", mbvFN)) ||
			
			!SetParam(bdiParamFN, "BDI", "USE_BDI_DATA", useBDIData) ||
			!SetParam(bdiParamFN, "BDI", "BDI_DATA", BDIDataFN) ||

			!SetParam(bdiParamFN, "BDI", "USE_ALT_OHA_DATA", useAltOHAData) ||
			!SetParam(bdiParamFN, "BDI", "ALT_OHA_DATA", altOHADataFN) ||

			!SetParam(bdiParamFN, "BDI", "SAR_Z", saRatio) ||
			!SetParam(bdiParamFN, "BDI", "SCALEHAB", scaleCond) ||
			!SetParam(bdiParamFN, "BDI", "MULTEHACOND", multEHACond) ||
			(useMBVAdd && !SetParam(bdiParamFN, "BDI", "MBV_ADD", mbvAdd)) ||
			(useMBVFac && !SetParam(bdiParamFN, "BDI", "MBV_FAC", mbvFac)))

			return msgErrorStr(-14, "parameters", bdiParamFN);
		
		std::cout << "Processing BDI parameter file: " << bdiParamFN << "\n";
		int result = CalcBDIForRegions(bdiParamFN);
		if (result != 0) return msgErrorStr(-15, "CalcBDIForRegions()", toStr(result));
	}

	//Run Benefits
	if (doBenefits) {
		for (auto benefit : benefits) {
			bool useAltGrids;
			std::string benFN, benMBVFN, altHab, altPrm, benDir;
			bool benMultFocal;
			float benFocalPower;
			float benSumPower;

			if (!GetParam(paramFN, benefit, "BEN_FN", benFN)) benFN = benefit;
			if (!GetParam(paramFN, benefit, "MBV_FN", benMBVFN)) benMBVFN = mbvFN;
			if (!GetFilesFolder(benFN, benDir)) benFN = workDir + pathSep + benFN;
			std::string benParamFN(benFN + ".par");

			std::cout << "Creating benefits parameter file: " << benParamFN << "\n";
			if (!GetParam(paramFN, benefit, "ALTGRIDS", useAltGrids)) return msgErrorStr(-1, "ALTGRIDS", paramFN);
			if (!GetParam(paramFN, benefit, "ALTHAB", altHab)) return msgErrorStr(-1, "ALTHAB", paramFN);
			if (!GetParam(paramFN, benefit, "ALTPRM", altPrm)) return msgErrorStr(-1, "ALTPRM", paramFN);

			if (!GetParam(paramFN, benefit, "MULTFOCAL", benMultFocal)) benMultFocal = multFocal;
			if (!GetParam(paramFN, benefit, "FOCALPOWER", benFocalPower)) benFocalPower = focalPower;
			if (!GetParam(paramFN, benefit, "SUMPOWER", benSumPower)) benSumPower = sumPower;

			if (!SetParam(benParamFN, "INPUTS", "HAB_FN", habFN)) return msgErrorStr(-14, "HAB_FN", benParamFN);
			if (!SetParam(benParamFN, "INPUTS", "PRM_FN", prmFN)) return msgErrorStr(-14, "PRM_FN", benParamFN);
			if (!SetParam(benParamFN, "INPUTS", "MBV_FN", benMBVFN)) return msgErrorStr(-14, "MBV_FN", benParamFN);
			if (!SetParam(benParamFN, "OUTPUTS", "BEN_FN", benFN)) return msgErrorStr(-14, "BEN_FN", benParamFN);

			if (!SetParam(benParamFN, "INPUTS", "ALTGRIDS", useAltGrids)) return msgErrorStr(-14, "ALTGRIDS", benParamFN);
			if (!SetParam(benParamFN, "INPUTS", "ALTHAB", altHab)) return msgErrorStr(-14, "ALTHAB", benParamFN);
			if (!SetParam(benParamFN, "INPUTS", "ALTPRM", altPrm)) return msgErrorStr(-14, "ALTPRM", benParamFN);

			if (!SetParam(benParamFN, "WINDOW", "LOOKUP_FN", lookupFN)) return msgErrorStr(-14, "LOOKUP_FN", benParamFN);
			if (!SetParam(benParamFN, "WINDOW", "TRANSLATE_FN", translateFN)) return msgErrorStr(-14, "TRANSLATE_FN", benParamFN);
			
			if (!SetParam(benParamFN, "ANALYSIS", "MULTFOCAL", benMultFocal)) return msgErrorStr(-14, "MULTFOCAL", benParamFN);
			if (!SetParam(benParamFN, "ANALYSIS", "FOCALPOWER", benFocalPower)) return msgErrorStr(-14, "FOCALPOWER", benParamFN);
			if (!SetParam(benParamFN, "ANALYSIS", "SUMPOWER", benSumPower)) return msgErrorStr(-14, "SUMPOWER", benParamFN);

			std::cout << "Processing benefits parameter file: " << benParamFN << "\n";
			int result = BenefitsCBA(benParamFN.c_str());
			if (result != 0) return msgErrorStr(-15, "BenefitsCBA()", toStr(result));
		}
	}
	return 0;
}