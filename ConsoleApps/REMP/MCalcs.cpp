/*
MCalcs.cpp - TODO Description
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


//MCalcs.cpp
//RunMCalcs() Creates MCalcs input grids then runs REMP CBA to calculate lambda and c
//TODO Allow geographic grids with parameters in metres

#include <string>
#include <cmath>
#include "MCalcs.h"
#include "FileUtilsLib.h"
#include "SpatialContextLib.h"

#ifdef _MSC_VER
#if _MSC_VER < 1910
#include <boost\filesystem.hpp>
namespace fileSys = boost::filesystem;
#elif _MSC_VER < 1920
#include <filesystem>
namespace fileSys = std::experimental::filesystem;
#else
#include <filesystem>
namespace fileSys = std::filesystem;
#endif
#else
#include <filesystem>
namespace fileSys = std::filesystem;
#endif

//Creates MCalcs input grids then runs REMP CBA to calculate lambda and c
int RunMCalcs(const char *paramFN)
{
	return RunMCalcs(std::string(paramFN));
}
int RunMCalcs(const std::string & paramFN)
{

	//Local variables
	std::string workDir, preFN;
	float lambda;
	double distMin, distMax, cellSize, mvh;
	bool rescaleMCalcsHab = true;

	//Get required parameters
	if (!GetParam(paramFN, "SPECIES", "MVH", mvh)) return msgErrorStr(-1, "MVH", paramFN);
	if (!GetParam(paramFN, "SPECIES", "CELLSIZE", cellSize)) return msgErrorStr(-1, "CELLSIZE", paramFN);
	if (!GetParam(paramFN, "SPECIES", "DISTMIN", distMin)) return msgErrorStr(-1, "DISTMIN", paramFN);
	if (!GetParam(paramFN, "SPECIES", "DISTMAX", distMax)) return msgErrorStr(-1, "DISTMAX", paramFN);
	if (!GetParam(paramFN, "ANALYSIS", "RESCALEHAB", rescaleMCalcsHab)) rescaleMCalcsHab = true;
	if (!GetParam(paramFN, "ANALYSIS", "PREFIX", preFN)) preFN = "";
	if (!GetParam(paramFN, "ANALYSIS", "WORKDIR", workDir) && !GetFilesFolder(paramFN, workDir)) workDir = ".";

	if (preFN.length() > 0 && preFN.back() != '_') preFN += '_';
	std::string habFN = workDir + "\\" + preFN + "MCalcs_Habitat";
	std::string prmFN = workDir + "\\" + preFN + "MCalcs_Permeability";
	std::string occFN = workDir + "\\" + preFN + "MCalcs_Occupancy";
	std::string mpcFN = workDir + "\\" + preFN + "MCalcs_MPC";

	//Create output folders if they don't already exist
	if (!(fileSys::is_directory(fileSys::path(workDir))) && !(fileSys::create_directories(fileSys::path(workDir))))
		return msgErrorStr(-99, "Unable to find or create folder:", workDir);

	//Calculte MCalcs input grid attributes
	//Approx cellSize in metres
	float cellArea = float(std::pow(cellSize, 2) / (cellSize > 1.0 ? 10000.0 : 0.000001));
	float cellSizeM = cellSize > 1.0 ? float(cellSize) : std::pow(cellArea, 0.5f);
	float prmMin = float(exp(-cellSizeM / distMin));
	float prmMax = float(exp(-cellSizeM / distMax));
	int grdDim = int(round(distMax * 64.0 / cellSizeM));
	grdDim += 1 - grdDim % 2;

	//Create MCalcs input grids
	if (!CreateCircleGrid(habFN, cellSize, mvh, grdDim, cellArea, 0.0f) ||
		!CreateCircleGrid(prmFN, cellSize, mvh, grdDim, prmMax, prmMin))
		return msgErrorStr(-99, "Unable to create MCALC input grids", habFN, prmFN);
	
	if (rescaleMCalcsHab) {
		double habSum = 0.0;
		if (!GetGridSum(habFN, habSum)) return msgErrorStr(-3, habFN);
		if (!ApplyLambdaToGrid(habFN, F2FLAMBDA(float(x * (mvh / habSum));))) return msgErrorStr(-4, habFN);
	}
	if (!SetParam(paramFN, "INPUTS", "HAB_FN0", habFN)) return msgErrorStr(-14, "HAB_FN0", paramFN);
	if (!SetParam(paramFN, "INPUTS", "PRM_FN0", prmFN)) return msgErrorStr(-14, "PRM_FN0", paramFN);
	if (!SetParam(paramFN, "OUTPUTS", "OCC_FN0", occFN)) return msgErrorStr(-14, "OCC_FN0", paramFN);
	if (!SetParam(paramFN, "OUTPUTS", "MPC_FN0", mpcFN)) return msgErrorStr(-14, "MPC_FN0", paramFN);

	//Run REMP CBA for MCalcs
	int result = OccupancyCBA_TS(paramFN.c_str());
	if (result != 0) return msgErrorStr(-15, "REMP_TS_CBA()", toStr(result));
	
	//Get lambda, c = e / Lambda
	if (!GetGridMax(mpcFN, lambda)) return msgErrorStr(-99, "Unable to get maximum value from", mpcFN);
	
	//Update values in MCalcs parameter file
	if (!SetParam(paramFN, "SPECIES", "LAMBDA", lambda))// ||
		return msgErrorStr(-14, "LAMBDA", paramFN);

	//Complete!
	msgText("RunMCalcs() Complete!");
	return 0;
}

/*
Create MVH grids at home range
Run NHA for HR MVH grids, gamma = max(NHA)
Resample HR MVH grids to dispersal range
Run MPC for DR MVH grids, lambda1 = MPC1 x gamma
Run Occ for DR MVH grids using labmda1, lambda2 = Si2 x gamma
Run Occ for DR MVH grids using labmda2, lambda3 = Si3 x gamma
*/

int RunMCalcsTriLambda(const std::string & paramFN)
{	
	//Local variables
	std::string hrLookupFN(""), hrTranslateFN(""), drLookupFN(""), drTranslateFN("");
	std::string workDir, preFN;

	float lambda, gamma;
	double hrDistMin, hrDistMax, hrCellSize;
	double drDistMin, drDistMax, drCellSize, mvh;
	int result;
	bool rescaleMCalcsHab = true;
	int hrIterations = 1;
	bool doHRIterations = true;
	bool useMaxPetalPerm = false;

	//Get required parameters
	if (!GetParam(paramFN, "WINDOW", "LOOKUP_FN", drLookupFN)) return msgErrorStr(-1, "LOOKUP_FN", paramFN);
	if (!GetParam(paramFN, "WINDOW", "TRANSLATE_FN", drTranslateFN)) return msgErrorStr(-1, "TRANSLATE_FN", paramFN);
	if (!GetParam(paramFN, "HRWINDOW", "LOOKUP_FN", hrLookupFN)) return msgErrorStr(-1, "LOOKUP_FN", paramFN);
	if (!GetParam(paramFN, "HRWINDOW", "TRANSLATE_FN", hrTranslateFN)) return msgErrorStr(-1, "TRANSLATE_FN", paramFN);
	if (!GetParam(paramFN, "SPECIES", "MVH", mvh)) return msgErrorStr(-1, "MVH", paramFN);
	if (!GetParam(paramFN, "SPECIES", "HRCELLSIZE", hrCellSize)) return msgErrorStr(-1, "HRCELLSIZE", paramFN);
	if (!GetParam(paramFN, "SPECIES", "HRDISTMIN", hrDistMin)) return msgErrorStr(-1, "HRDISTMIN", paramFN);
	if (!GetParam(paramFN, "SPECIES", "HRDISTMAX", hrDistMax)) return msgErrorStr(-1, "HRDISTMAX", paramFN);
	if (!GetParam(paramFN, "SPECIES", "CELLSIZE", drCellSize)) return msgErrorStr(-1, "CELLSIZE", paramFN);
	if (!GetParam(paramFN, "SPECIES", "DISTMIN", drDistMin)) return msgErrorStr(-1, "DISTMIN", paramFN);
	if (!GetParam(paramFN, "SPECIES", "DISTMAX", drDistMax)) return msgErrorStr(-1, "DISTMAX", paramFN);
	if (!GetParam(paramFN, "ANALYSIS", "RESCALEHAB", rescaleMCalcsHab)) rescaleMCalcsHab = true;
	if (!GetParam(paramFN, "ANALYSIS", "PREFIX", preFN)) preFN = "";
	if (!GetParam(paramFN, "ANALYSIS", "WORKDIR", workDir) && !GetFilesFolder(paramFN, workDir)) workDir = ".";

	//if (!GetParam(paramFN, "ANALYSIS", "DOHRITERATIONS", doHRIterations)) doHRIterations = false;
	if (!GetParam(paramFN, "ANALYSIS", "HRITERATIONS", hrIterations)) doHRIterations = false;

	if (!GetParam(paramFN, "ANALYSIS", "MAXPETALPERM", useMaxPetalPerm)) useMaxPetalPerm = false;

	
	if (preFN.length() > 0 && preFN.back() != '_') preFN += '_';
	std::string hrDir = workDir + "\\Homerange";
	std::string drDir = workDir + "\\Dispersal";
	std::string hrHabFN = hrDir + "\\" + preFN + "MCalcs_HR_Habitat";
	std::string hrPrmFN = hrDir + "\\" + preFN + "MCalcs_HR_Permeability";
	std::string hrNHAFN = hrDir + "\\" + preFN + "MCalcs_HR_NHA";
	std::string hrMPCFN = hrDir + "\\" + preFN + "MCalcs_HR_MPC";
	std::string drHabFN = drDir + "\\" + preFN + "MCalcs_DR_Habitat";
	std::string drPrmFN = drDir + "\\" + preFN + "MCalcs_DR_Permeability";
	std::string drNHAFN = drDir + "\\" + preFN + "MCalcs_DR_NHA";
	std::string drOccFN = drDir + "\\" + preFN + "MCalcs_DR_Occupancy";
	std::string drMPCFN = drDir + "\\" + preFN + "MCalcs_DR_MPC";

	//Create output folders if they don't already exist
	if (!(fileSys::is_directory(fileSys::path(workDir))) && !(fileSys::create_directories(fileSys::path(workDir))))
		return msgErrorStr(-99, "Unable to find or create folder:", workDir);
	if (!(fileSys::is_directory(fileSys::path(hrDir))) && !(fileSys::create_directories(fileSys::path(hrDir))))
		return msgErrorStr(-99, "Unable to find or create folder:", hrDir);
	if (!(fileSys::is_directory(fileSys::path(drDir))) && !(fileSys::create_directories(fileSys::path(drDir))))
		return msgErrorStr(-99, "Unable to find or create folder:", drDir);
	
	//TODO cellSizeInMetres

	//Get cell area and permeability values for hr and dr scales
	float hrCellArea = float(std::pow(hrCellSize, 2) / (hrCellSize > 1.0 ? 10000.0 : 0.000001));
	float drCellArea = float(std::pow(drCellSize, 2) / (drCellSize > 1.0 ? 10000.0 : 0.000001));
	float hrPrmMin = float(std::exp(-hrCellSize / hrDistMin));
	float hrPrmMax = float(std::exp(-hrCellSize / hrDistMax));
	//float drPrmMin = float(exp(-drCellSize / drDistMin));
	//float drPrmMax = float(exp(-drCellSize / drDistMax));
	
	int hrGrdDim = int(round(hrDistMax * 64.0 / hrCellSize) * (drCellSize / hrCellSize));
	int drGrdDim = int(round(drDistMax * 64.0 / drCellSize));
	hrGrdDim += 1 - hrGrdDim % 2;
	drGrdDim += 1 - drGrdDim % 2;
	float tmpMin = 0.0f, tmpMax = 0.0f;

	//Create HR MVH habitat and permeability grids at home range
	if (!CreateCircleGrid(hrHabFN, hrCellSize, mvh, hrGrdDim, hrCellArea, 0.0f) ||
		!CreateCircleGrid(hrPrmFN, hrCellSize, mvh, hrGrdDim, hrPrmMax, hrPrmMin))
		return msgErrorStr(-99, "Unable to create MCALC input grids\n", hrHabFN, hrPrmFN);
	
	//Rescale HR MVH habitat so that cell values sum exactly to MVH
	if (rescaleMCalcsHab) {
		double habSum = 0.0;
		if (!GetGridSum(hrHabFN, habSum) || habSum == 0.0) return msgErrorStr(-3, hrHabFN);
		if (!ApplyLambdaToGrid(hrHabFN, F2FLAMBDA(float(x * (mvh / habSum));))) return msgErrorStr(-4, hrHabFN);
	}

	std::string hrParamFN(hrDir + "\\" + preFN + "homerange.par");
	if (doHRIterations) {
		if (!SetParam(hrParamFN, "INPUTS", "HAB_FN0", hrHabFN) ||
			!SetParam(hrParamFN, "INPUTS", "PRM_FN0", hrPrmFN) ||
			!SetParam(hrParamFN, "OUTPUTS", "MPC_FN0", hrMPCFN) ||
			!SetParam(hrParamFN, "WINDOW", "LOOKUP_FN", hrLookupFN) ||
			!SetParam(hrParamFN, "WINDOW", "TRANSLATE_FN", hrTranslateFN) ||
			!SetParam(hrParamFN, "SPECIES", "LAMBDA", 1.0) || //REMP_TS_CBA() wants LAMBDA but it's not needed here, make optional
			!SetParam(hrParamFN, "SPECIES", "#NOTE", "LAMBDA is not used here but parameter needed by REMP CBA") ||
			!SetParam(hrParamFN, "ANALYSIS", "ITERATIONS", hrIterations) ||
			!SetParam(hrParamFN, "ANALYSIS", "MAXPETALPERM", useMaxPetalPerm))
			return msgErrorStr(-14, "parameters", hrParamFN);

		std::cout << "Performing homerange MPC analysis using parameters from:\n" << hrParamFN << "\n\n";
		result = OccupancyCBA_TS(hrParamFN.c_str());
		if (result != 0) return msgErrorStr(-15, "REMP_TS_CBA()", toStr(result));
		hrNHAFN = hrMPCFN;
	}
	else {
		//Set HR parameters then run NHA for HR MVH grids to get gamma = max(NHA)
		if (!SetParam(hrParamFN, "INPUTS", "HAB_FN", hrHabFN) ||
			!SetParam(hrParamFN, "INPUTS", "PRM_FN", hrPrmFN) ||
			!SetParam(hrParamFN, "OUTPUTS", "NHA_FN", hrNHAFN) ||
			!SetParam(hrParamFN, "WINDOW", "LOOKUP_FN", hrLookupFN) ||
			!SetParam(hrParamFN, "WINDOW", "TRANSLATE_FN", hrTranslateFN))
			return msgErrorStr(-14, "parameters", hrParamFN);

		result = ContextCBA(hrParamFN.c_str());
		if (result != 0) return msgErrorStr(-15, "ContextCBA()", toStr(result));
	}
	if (!GetGridMax(hrNHAFN, gamma)) return msgErrorStr(-3, hrNHAFN);

	//Resample HR MVH habitat to dispersal range and rescale
	if (!ResampleGrid(hrHabFN, drHabFN, drCellSize, 5)) return msgErrorStr(-13, drHabFN);
	if (!ApplyLambdaToGrid(drHabFN, F2FLAMBDA(float(drCellArea) * x / hrCellArea;))) return msgErrorStr(-13, drHabFN);
	
	//Rescale DR MVH habitat so that cell values sum exactly to MVH
	if (rescaleMCalcsHab) {
		double habSum = 0.0;
		if (!GetGridSum(drHabFN, habSum) || habSum == 0.0) return msgErrorStr(-3, drHabFN);
		if (!ApplyLambdaToGrid(drHabFN, F2FLAMBDA(float(x * (mvh / habSum));))) return msgErrorStr(-4, drHabFN);
	}
	
	//Resample HR MVH permeability to dispersal range and rescale
	if (!CreatePrmGrid(drHabFN, drPrmFN, drDistMin, drDistMax)) return msgErrorStr(-13, drPrmFN);
	//if (!ResampleGrid(hrPrmFN, drPrmFN, drHabFN, 5)) return msgErrorStr(-13, drPrmFN);
	//if (!GetGridMinMax(drPrmFN, tmpMin, tmpMax)) return msgErrorStr(-3, drPrmFN);
	//if (!ApplyLambdaToGrid(drPrmFN, F2FLAMBDA(((x - tmpMin) / (tmpMax - tmpMin)) * float(drPrmMax - drPrmMin) + float(drPrmMin);))) return msgErrorStr(-4, drPrmFN);

	//Resample HR MVH NHA to dispersal range
	if (!ResampleGrid(hrNHAFN, drNHAFN, drHabFN, 5)) return msgErrorStr(-13, drNHAFN);

	//Set DR input filename parameters
	if (!SetParam(paramFN, "INPUTS", "HAB_FN0", drHabFN)) return msgErrorStr(-14, "HAB_FN0", paramFN);
	if (!SetParam(paramFN, "INPUTS", "PRM_FN0", drPrmFN)) return msgErrorStr(-14, "PRM_FN0", paramFN);
	if (!SetParam(paramFN, "INPUTS", "NHA_FN0", drNHAFN)) return msgErrorStr(-14, "NHA_FN0", paramFN);

	//Run MPC for DR MVH grids, lambda1 = max(MPC) * G
	//Run Occ for DR MVH grids using labmda1, lambda2 = max(Si) * G
	//Run Occ for DR MVH grids using labmda2, lambda3 = max(Si) * G
	for (int s = 1; s < 4; s++) {

		if (!SetParam(paramFN, "OUTPUTS", "OCC_FN0", drOccFN + toStr(s))) return msgErrorStr(-14, "OCC_FN0", paramFN);
		if (!SetParam(paramFN, "OUTPUTS", "MPC_FN0", drMPCFN + toStr(s))) return msgErrorStr(-14, "MPC_FN0", paramFN);

		result = OccupancyCBA_TS(paramFN.c_str());
		if (result != 0) return msgErrorStr(-15, "REMP_TS_CBA()", toStr(result));

		if (s == 1) { if (!GetGridMax(drMPCFN + toStr(s), lambda)) return msgErrorStr(-3, drMPCFN + toStr(s)); }
		else if (!GetGridMax(drOccFN + toStr(s) + "_Si", lambda)) return msgErrorStr(-3, drOccFN + toStr(s) + "_Si");
		lambda *= gamma;

		if (!SetParam(paramFN, "SPECIES", "LAMBDA", lambda) ||
			!SetParam(paramFN, "SPECIES", "LAMBDA" + toStr(s), lambda))
			return msgErrorStr(-14, "LAMBDA and LAMBDA" + toStr(s), paramFN);
	}

	//Set lambda to lambda2 for region grid occupancy
	if (!GetParam(paramFN, "SPECIES", "LAMBDA2", lambda) ||
		!SetParam(paramFN, "SPECIES", "LAMBDA", lambda))
		return msgErrorStr(-14, "LAMBDA and C", paramFN);

	//Complete!
	msgText("RunMCalcsTriLambda() Complete!");
	return 0;
}

/**/