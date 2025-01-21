/*
REMP.cpp - TODO Description
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


#include <fstream>
#include <iostream>
#include <cmath>
#include "FileUtilsLib.h"
#include "SpatialContextLib.h"
#include "MCalcs.h"

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


//TODO Allow geographic grids with parameters in metres
//TODO Rename grids to match equations
//TODO FIX using tif causes bug when no data value differs between habitat and permeability grids.

/*
Includes new triLambda MCalcs version
Create MVH grids at home range
Run NHA for HR MVH grids, gamma = max(NHA)
Resample HR MVH grids to dispersal range
Run MPC for DR MVH grids, lambda1 = MPC1 x gamma
Run Occ for DR MVH grids using labmda1, lambda2 = Si2 x gamma
Run Occ for DR MVH grids using labmda2, lambda3 = Si3 x gamma
Resample region grids to HR if HR!=Source
Run NHA for HR region grids
Resample region grids to DR
Run Occ for DR region grids using lambda2, Pu = Si / (Si + lambda3 / NHA)
*/

//Disable sprintf is depreciated
#pragma warning(disable : 4996)

using uint = unsigned int;

int RunREMP(const std::string & srcParamFN);

int main(int argc, char *argv[])
{
	int result = 0;
	//Options
	if (argc != 2) {
		std::cout << "\nUsage: " << argv[0] << " ParameterFile.par \n"
			<< "\nIf ParameterFile.par exists it is processed, otherwise it is created as a template.\n\n";
		return -1;
	}
	std::string argv1(argv[1]);
	if (!FileExists(argv1)) {
		std::cout << "\nThe parameter file " << argv1 << " does not exist and is being created as a template\n\n";
		//Input parameters
		if (!SetParam(argv[1], "INPUTS", "HAB_FN0", "Filename. Habitat input grid file name without extension for time step 0. Additional _FN1..._FNt-1 required for t time steps.") ||
		    !SetParam(argv[1], "INPUTS", "PRM_FN0", "Filename. Permeability input grid file name without extension for time step 0. Additional _FN1..._FNt-1 required for t time steps. Optional, default derived from HAB_FNt.") ||
		    !SetParam(argv[1], "INPUTS", "NHA_FN0", "Filename. Existing source or home range NHA file name if supplying. Additional _FN1..._FNt-1 required for t time steps. Optional, default created at home range.") ||
			!SetParam(argv[1], "INPUTS", "SCALEHABHA", "true/false, Whether to scale habitat input to area in hectares. Optional, default = true.") ||
			!SetParam(argv[1], "INPUTS", "SRCHABMAX", "Floating point number. Maximum possible source habitat value used to scale the habitat input. Optional, default derived from HAB_FNt.") ||

		//Window parameters
			!SetParam(argv[1], "WINDOW", "#NOTE", "Use [WINDOW] section if home and dispersal range window parameters match, else use seperate [HRWINDOW] and [DRWINDOW] sections.") ||
			!SetParam(argv[1], "WINDOW", "SRCSIZE", "Integer. The source size of the home and dispersal range search window in cells. Optional, default use existing lookup and translate tables.") ||
			!SetParam(argv[1], "WINDOW", "DSTSIZE", "Integer. The destination size of the home and dispersal range search window in petals. Optional, default use existing lookup and translate tables.") ||
			!SetParam(argv[1], "WINDOW", "ZONERATIO", "Floating point number. The area ratio increase of each petal outwards. Must be >= 1.0. Optional, default use existing lookup and translate tables.") ||
			!SetParam(argv[1], "WINDOW", "SAVETMPFILES", "true/false. whether to save temporary petal files when generating lookup and translate tables. Optional, default = true.") ||
			!SetParam(argv[1], "WINDOW", "LOOKUP_FN", "Filename. The petals lookup table file name with extension. Optional, default created using src/dst size and zone ratio.") ||
			!SetParam(argv[1], "WINDOW", "TRANSLATE_FN", "Filename. The petals translate table filename with extension. Optional, default created using src/dst size and zone ratio.") ||

		//Species parameters
		    !SetParam(argv[1], "SPECIES", "MVH", "Floating point number. Minimum viable habitat in hectares.") ||
		    !SetParam(argv[1], "SPECIES", "HRCELLSIZE", "Floating point number. Home range analysis cell size.") ||
		    !SetParam(argv[1], "SPECIES", "DRCELLSIZE", "Floating point number. Dispersal range analysis cell size.") ||
		    !SetParam(argv[1], "SPECIES", "HRDISTMIN", "Floating point number. Minimum home range movement distance") ||
		    !SetParam(argv[1], "SPECIES", "HRDISTMAX", "Floating point number. Maximum home range movement distance") ||
		    !SetParam(argv[1], "SPECIES", "DRDISTMIN", "Floating point number. Minimum dispersal range movement distance") ||
		    !SetParam(argv[1], "SPECIES", "DRDISTMAX", "Floating point number. Maximum dispersal range movement distance") ||
		    !SetParam(argv[1], "SPECIES", "LAMBDA", "Floating point number. Lambda parameter or optionally LAMBDA1 LAMBDA2 LAMBDA3 parameters. Optional, default derived from MCalcs.") ||
		    !SetParam(argv[1], "SPECIES", "C", "Floating point number. DEPRECIATED Use LAMBDA! Colonisation parameter c. Optional, default derived from MCalcs.") ||
		    !SetParam(argv[1], "SPECIES", "E", "Floating point number. DEPRECIATED Use LAMBDA! Extinction parameter e. Optional, default = 1.0.") ||

		//Analysis parameters
			!SetParam(argv[1], "ANALYSIS", "WORKDIR", "Folder name. Path to an existing folder where outputs are written to. Optional, default = parameter file's folder if resolvable.") ||
			!SetParam(argv[1], "ANALYSIS", "PREFIX", "Filename. Prefix that's added to the beginning of each generated file. Optional, default = not used.") ||
			!SetParam(argv[1], "ANALYSIS", "TRILAMBDA", "true/false. Whether to generate three lambdas during MCalcs using lambda2 for occ and lambda3 for Pu. Optional, default = true.") ||
			!SetParam(argv[1], "ANALYSIS", "RESCALEHAB", "true/false. Whether to rescale MCalcs habitat values so they sum to MVH more precisely. Optional, default = true.") ||
		    !SetParam(argv[1], "ANALYSIS", "FINALPIATHR", "true/false. Whether to calculate Pi at home range rather than source resolution. Optional, default = false.") || 
		    !SetParam(argv[1], "ANALYSIS", "TIMESTEPS", "Integer. Number of time steps to be processed. Must provide _FN1...FNt-1 Filenames. Optional, default = 1.") ||
			!SetParam(argv[1], "ANALYSIS", "ITERATIONS", "Integer. Number of CBA iterations to perform (per time step). Optional, default = 1.") ||
			!SetParam(argv[1], "ANALYSIS", "HRITERATIONS", "Integer. Number of home range CBA iterations to perform (per time step) for home range MPC. Optional, default use home range NHA with no iterations.") ||
			
			!SetParam(argv[1], "ANALYSIS", "COUPLED", "true/false. Whether occupancy is coupled between time steps. Optional, default = false, need OCC_FN when true.") ||
			!SetParam(argv[1], "ANALYSIS", "DEBUG", "true/false. Whether intermediate files are written each CBA iteration for small datasets (always written for large datasets). Optional, default = false.") ||

			!SetParam(argv[1], "ANALYSIS", "MAXPETALPERM", "true/false. Whether to use the maximum permeability in each petal rather than the average, default = false.") ||

			!SetParam(argv[1], "ANALYSIS", "DOMCALCS", "true/false. Whether MCalcs analysis is performed. Optional, default = true.") ||
			!SetParam(argv[1], "ANALYSIS", "ONLYMCALCS","true/false. Whether only MCalcs analysis is performed. Optional, default = false.") ||
			!SetParam(argv[1], "ANALYSIS", "DOSRINPUTS", "true/false. Whether source resolution grids are generated. Optional, default = true.") ||
			!SetParam(argv[1], "ANALYSIS", "DOHRINPUTS", "true/false. Whether homerange resolution inputs are generated. Optional, default = true.") ||
			!SetParam(argv[1], "ANALYSIS", "DODRINPUTS", "true/false. Whether dispersal resolution inputs are generated. Optional, default = true.") ||
			!SetParam(argv[1], "ANALYSIS", "DOHOMERANGE", "true/false. Whether homerange analysis is performed. Optional, default = true.") ||
			!SetParam(argv[1], "ANALYSIS", "DODISPERSAL", "true/false. Whether homerange analysis is performed. Optional, default = true.") ||
			!SetParam(argv[1], "ANALYSIS", "DOFINALPI", "true/false. Whether final occupancy analysis is performed. Optional, default = true.") ||
			
			!SetParam(argv[1], "ANALYSIS", "CLEANALL", "true/false. Whether All intermediate folders and files are removed once processed. Optional, default = false.") ||
			!SetParam(argv[1], "ANALYSIS", "CLEANMCALCS", "true/false. Whether MCalcs folder and files are removed once processed. Optional, default = false.") ||
			!SetParam(argv[1], "ANALYSIS", "CLEANSR", "true/false. Whether Source folder and files are removed once processed. Optional, default = false.") ||
			!SetParam(argv[1], "ANALYSIS", "CLEANHR", "true/false. Whether Homerange folder and files are removed once processed. Optional, default = false.") ||
			!SetParam(argv[1], "ANALYSIS", "CLEANDR", "true/false. Whether Dispersal folder and files are removed once processed. Optional, default = false."))

			return msgErrorStr(-99, "Unable to create parameter file template:", argv1);
		return -2;
	}

	result = RunREMP(argv1);
	if (result != 0) std::cout << "Error code: " << result << "\n";
	else std::cout << "Complete!" << "\n";
	return result;
}



int RunREMP(const std::string & srcParamFN) 
{
	int result = 0;
	std::string srcParamDir(""), workDir("");
	if (!FileExists(srcParamFN)) return msgErrorStr(-99, "Can't find parameter file:", srcParamFN);
	if (!GetFilesFolder(srcParamFN, srcParamDir)) srcParamDir = ".";

	////////////////////////////////////////////////////////////////////////////////////
	//Get parameters from source parameter file
	////////////////////////////////////////////////////////////////////////////////////
	//Input parameters
	std::string srHabFN0(""), preFN = "";
	bool applyAreaToHab = true;
	bool haveSrcHabMax = true;
	float srHabMax = 0.0f;

	//Species parameters
	double mvh, srCellSize, hrCellSize, drCellSize;
	float hrDistMin, hrDistMax, drDistMin, drDistMax, lambda, lambda1, lambda2, lambda3;

	//Analysis parameters
	bool doTriLambda = false;
	bool rescaleMCalcsHab = true;
	bool doFinalStepAtHR = false;
	int timeSteps = 1;
	int iterations = 1;

	int hrIterations = 1;
	bool doHRIterations = true;

	bool coupled = false;
	bool rempDebug = false;

	bool doMCalcs = true;
	bool onlyMCalcs = false;
	bool doSRInputs = true;
	bool doHRInputs = true;
	bool doDRInputs = true;
	bool doHomerange = true;
	bool doDispersal = true;
	bool doFinalPi = true;

	bool useMaxPetalPerm = false;

	bool cleanAll = false;
	bool cleanMC = false;
	bool cleanSR = false;
	bool cleanHR = false;
	bool cleanDR = false;

	//Input parameters
	if (!GetParam(srcParamFN, "INPUTS", "HAB_FN0", srHabFN0)) return msgErrorStr(-1, "HAB_FN0", srcParamFN);
	if (!GetParam(srcParamFN, "INPUTS", "SCALEHABHA", applyAreaToHab)) applyAreaToHab = true;
	if (!GetParam(srcParamFN, "INPUTS", "SRCHABMAX", srHabMax)) haveSrcHabMax = false;

	//Species parameters
	if (!GetParam(srcParamFN, "SPECIES", "MVH", mvh)) return msgErrorStr(-1, "MVH", srcParamFN);
	if (!GetParam(srcParamFN, "SPECIES", "HRCELLSIZE", hrCellSize)) return msgErrorStr(-1, "HRCELLSIZE", srcParamFN);
	if (!GetParam(srcParamFN, "SPECIES", "DRCELLSIZE", drCellSize)) return msgErrorStr(-1, "DRCELLSIZE", srcParamFN);
	if (!GetParam(srcParamFN, "SPECIES", "HRDISTMIN", hrDistMin)) return msgErrorStr(-1, "HRDISTMIN", srcParamFN);
	if (!GetParam(srcParamFN, "SPECIES", "HRDISTMAX", hrDistMax)) return msgErrorStr(-1, "HRDISTMAX", srcParamFN);
	if (!GetParam(srcParamFN, "SPECIES", "DRDISTMIN", drDistMin)) return msgErrorStr(-1, "DRDISTMIN", srcParamFN);
	if (!GetParam(srcParamFN, "SPECIES", "DRDISTMAX", drDistMax)) return msgErrorStr(-1, "DRDISTMAX", srcParamFN);
	if (!GetParam(srcParamFN, "SPECIES", "LAMBDA", lambda)) lambda = 1.0f;
	if (!GetParam(srcParamFN, "SPECIES", "LAMBDA1", lambda1)) lambda1 = 1.0f;
	if (!GetParam(srcParamFN, "SPECIES", "LAMBDA2", lambda2)) lambda2 = 1.0f;
	if (!GetParam(srcParamFN, "SPECIES", "LAMBDA3", lambda3)) lambda3 = 1.0f;

	//Analysis parameters
	if (!GetParam(srcParamFN, "ANALYSIS", "WORKDIR", workDir)) workDir = srcParamDir;
	if (!GetParam(srcParamFN, "ANALYSIS", "PREFIX", preFN)) preFN = "";
	if (!GetParam(srcParamFN, "ANALYSIS", "TRILAMBDA", doTriLambda)) doTriLambda = true;
	if (!GetParam(srcParamFN, "ANALYSIS", "RESCALEHAB", rescaleMCalcsHab)) rescaleMCalcsHab = true;
	if (!GetParam(srcParamFN, "ANALYSIS", "FINALPIATHR", doFinalStepAtHR)) doFinalStepAtHR = false;
	if (!GetParam(srcParamFN, "ANALYSIS", "TIMESTEPS", timeSteps)) timeSteps = 1;
	if (!GetParam(srcParamFN, "ANALYSIS", "ITERATIONS", iterations)) iterations = 1;

	if (!GetParam(srcParamFN, "ANALYSIS", "HRITERATIONS", hrIterations)) doHRIterations = false;

	if (!GetParam(srcParamFN, "ANALYSIS", "MAXPETALPERM", useMaxPetalPerm)) useMaxPetalPerm = false;

	if (!GetParam(srcParamFN, "ANALYSIS", "COUPLED", coupled)) coupled = false;
	if (!GetParam(srcParamFN, "ANALYSIS", "DEBUG", rempDebug)) rempDebug = false;

	if (!GetParam(srcParamFN, "ANALYSIS", "DOMCALCS", doMCalcs)) doMCalcs = true;
	if (!GetParam(srcParamFN, "ANALYSIS", "ONLYMCALCS", onlyMCalcs)) onlyMCalcs = false;
	if (!GetParam(srcParamFN, "ANALYSIS", "DOSRINPUTS", doSRInputs)) doSRInputs = true;
	if (!GetParam(srcParamFN, "ANALYSIS", "DOHRINPUTS", doHRInputs)) doHRInputs = true;
	if (!GetParam(srcParamFN, "ANALYSIS", "DODRINPUTS", doDRInputs)) doDRInputs = true;
	if (!GetParam(srcParamFN, "ANALYSIS", "DOHOMERANGE", doHomerange)) doHomerange = true;
	if (!GetParam(srcParamFN, "ANALYSIS", "DODISPERSAL", doDispersal)) doDispersal = true;
	if (!GetParam(srcParamFN, "ANALYSIS", "DOFINALPI", doFinalPi)) doFinalPi = true;

	if (!GetParam(srcParamFN, "ANALYSIS", "CLEANALL", cleanAll)) cleanAll = false;
	if (!GetParam(srcParamFN, "ANALYSIS", "CLEANMCALCS", cleanMC)) cleanMC = false;
	if (!GetParam(srcParamFN, "ANALYSIS", "CLEANSR", cleanSR)) cleanSR = false;
	if (!GetParam(srcParamFN, "ANALYSIS", "CLEANHR", cleanHR)) cleanHR = false;
	if (!GetParam(srcParamFN, "ANALYSIS", "CLEANDR", cleanDR)) cleanDR = false;

	//Output folders and parameter filenames
	if (preFN.length() > 0 && preFN.back() != '_') preFN += '_';
	std::string srDir = workDir + "\\Source";
	std::string mcDir = workDir + "\\MCalcs";
	std::string hrDir = workDir + "\\Homerange";
	std::string drDir = workDir + "\\Dispersal";
	std::string ocDir = workDir + "\\Occupancy";
	std::string mcParamFN = workDir + "\\" + preFN + "mcalcs.par";
	std::string hrParamFN = workDir + "\\" + preFN + "homerange.par";
	std::string drParamFN = workDir + "\\" + preFN + "dispersal.par";
	std::string tmpFN = workDir + "\\tempData";

	//Create output folders if they don't already exist
	std::cout << "Creating work directory and folder structure\n";
	if (!(fileSys::is_directory(fileSys::path(workDir))) && !(fileSys::create_directories(fileSys::path(workDir))))
		return msgErrorStr(-99, "Unable to find or create folder:", workDir);
	if (!(fileSys::is_directory(fileSys::path(mcDir))) && !(fileSys::create_directories(fileSys::path(mcDir))))
		return msgErrorStr(-99, "Unable to find or create folder:", mcDir);
	if (!(fileSys::is_directory(fileSys::path(srDir))) && !(fileSys::create_directories(fileSys::path(srDir))))
		return msgErrorStr(-99, "Unable to find or create folder:", srDir);
	if (!(fileSys::is_directory(fileSys::path(hrDir))) && !(fileSys::create_directories(fileSys::path(hrDir))))
		return msgErrorStr(-99, "Unable to find or create folder:", hrDir);
	if (!(fileSys::is_directory(fileSys::path(drDir))) && !(fileSys::create_directories(fileSys::path(drDir))))
		return msgErrorStr(-99, "Unable to find or create folder:", drDir);
	if (!(fileSys::is_directory(fileSys::path(ocDir))) && !(fileSys::create_directories(fileSys::path(ocDir))))
		return msgErrorStr(-99, "Unable to find or create folder:", ocDir);

	//Read in tif habitat grids. Hack job for SoS PLP project
	//Note using tif causes bug when no data value differs between habitat and permeability grids.
	if (GetFileExt(srHabFN0).compare("tif") == 0) {
		igstream inTif;
		ogstream outFlt;

		std::cout << "Copying input tif habitat grid to floating point grid\n";
		inTif.open(srHabFN0);
		if (!inTif.is_open()) return msgErrorStr(-3, srHabFN0);
		if (inTif.dataType() != 6) return msgErrorStr(-99, "Input tif habitat grid data type is not float32");

		outFlt.open(GridFN(srDir + "\\" + preFN + "SR_Habitat_t0"));
		outFlt.copyHeader(inTif);
		if (!outFlt.is_open()) return msgErrorStr(-3, srHabFN0);

		std::unique_ptr<float[]> data = std::make_unique<float[]>(inTif.nCols());

		for (uint r = 0; r < inTif.nRows(); r++) {
			inTif.read((char*)data.get(), sizeof(float) * inTif.nCols());
			outFlt.write((const char*)data.get(), sizeof(float) * inTif.nCols());
		}
		inTif.close();
		outFlt.close();
		srHabFN0 = (srDir + "\\" + preFN + "SR_Habitat_t0");
	}
	//End of hack



	//Check habitat max, cellsize and distances
	std::cout << "Getting maximum source resolution habitat value: ";
	if (!haveSrcHabMax && !GetGridMaxPow10(srHabFN0, srHabMax) && srHabMax < 1.0f) return msgErrorStr(-99, "Unable to get source habitat max from", srHabFN0);
	std::cout << srHabMax << "\n\n";

	if (!GetHeader(srHabFN0, "cellsize", srCellSize)) return msgErrorStr(-5, srHabFN0);
	if (hrCellSize < srCellSize) return msgErrorStr(-99, "Home range cell size must be greater or equal to source cell size");
	if (hrDistMin > hrDistMax) return msgErrorStr(-99, "Minimum home range distance is greater than maximum");
	if (drDistMin > drDistMax) return msgErrorStr(-99, "Minimum dispersal distance is greater than maximum");
	if (hrDistMin == 0.0f || hrDistMax == 0.0f || drDistMin == 0.0f || drDistMax == 0.0f)
		return msgErrorStr(-99, "Homerange and dispersal distances must be greater than 0");

	//Calculate cell area values for src, hr and dr scales
	//TODO cellSizeInMetres
	//TODO better to use projection than cellsize when switching to gstreams
	double srCellArea = std::pow(srCellSize, 2) / (srCellSize > 1.0f ? 10000.0 : 0.000001);
	double hrCellArea = std::pow(hrCellSize, 2) / (hrCellSize > 1.0f ? 10000.0 : 0.000001);
	double drCellArea = std::pow(drCellSize, 2) / (drCellSize > 1.0f ? 10000.0 : 0.000001);

	//Get petal lookup and translate file names or create them if window parameters provided
	std::string hrLookupFN(""), hrTranslateFN(""), drLookupFN(""), drTranslateFN("");
	std::cout << "Creating petal lookup and translate tables\n";
	if (ParamSectionExists(srcParamFN, "WINDOW") && GetPetalFiles(srcParamFN, "WINDOW", workDir, hrLookupFN, hrTranslateFN) == 0) {
		drLookupFN = hrLookupFN;
		drTranslateFN = hrTranslateFN;
	}
	else if (!ParamSectionExists(srcParamFN, "HRWINDOW") ||
		!ParamSectionExists(srcParamFN, "DRWINDOW") ||
		GetPetalFiles(srcParamFN, "HRWINDOW", workDir, hrLookupFN, hrTranslateFN) != 0 ||
		GetPetalFiles(srcParamFN, "DRWINDOW", workDir, drLookupFN, drTranslateFN) != 0)
		return msgErrorStr(-99, "Unable to get home and dispersal range petals from", srcParamFN);

	//Run MCalcs using MVH and Min/Max habitat and 1/a values to calculate lambda
	if (doMCalcs) {
		msgText("Running MCalcs");

		//Create the MCalcs parameter file
		if (!SetParam(mcParamFN, "WINDOW", "LOOKUP_FN", drLookupFN) ||
			!SetParam(mcParamFN, "WINDOW", "TRANSLATE_FN", drTranslateFN) ||
			!SetParam(mcParamFN, "SPECIES", "MVH", mvh) ||
			!SetParam(mcParamFN, "SPECIES", "CELLSIZE", drCellSize) ||
			!SetParam(mcParamFN, "SPECIES", "DISTMIN", drDistMin) ||
			!SetParam(mcParamFN, "SPECIES", "DISTMAX", drDistMax) ||
			!SetParam(mcParamFN, "SPECIES", "LAMBDA", 1.0f) ||
			!SetParam(mcParamFN, "ANALYSIS", "WORKDIR", mcDir) ||
			!SetParam(mcParamFN, "ANALYSIS", "PREFIX", preFN) ||
			!SetParam(mcParamFN, "ANALYSIS", "RESCALEHAB", rescaleMCalcsHab) ||
			!SetParam(mcParamFN, "ANALYSIS", "TIMESTEPS", 1) ||
			!SetParam(mcParamFN, "ANALYSIS", "ITERATIONS", iterations) ||
			!SetParam(mcParamFN, "ANALYSIS", "COUPLED", coupled) ||
			!SetParam(mcParamFN, "ANALYSIS", "MAXPETALPERM", useMaxPetalPerm) ||
			!SetParam(mcParamFN, "ANALYSIS", "DEBUG", rempDebug))
			return msgErrorStr(-14, "parameters", mcParamFN);

		if (doTriLambda) {
			if (!SetParam(mcParamFN, "HRWINDOW", "LOOKUP_FN", hrLookupFN) ||
				!SetParam(mcParamFN, "HRWINDOW", "TRANSLATE_FN", hrTranslateFN) ||
				!SetParam(mcParamFN, "SPECIES", "HRCELLSIZE", hrCellSize) ||
				!SetParam(mcParamFN, "SPECIES", "HRDISTMIN", hrDistMin) ||
				!SetParam(mcParamFN, "SPECIES", "HRDISTMAX", hrDistMax))
				return msgErrorStr(-14, "parameters", mcParamFN);
			
			if (doHRIterations && !SetParam(mcParamFN, "ANALYSIS", "HRITERATIONS", hrIterations))
				return msgErrorStr(-14, "parameters", mcParamFN);

			//Run MCalcs generating three lambdas
			result = RunMCalcsTriLambda(mcParamFN);
			if (result != 0) return msgErrorStr(-15, "RunMCalcsTriLambda()", toStr(result));

			//Get lamda parameters from MCalcs parameter file and set in source parameter file
			if (!GetParam(mcParamFN, "SPECIES", "LAMBDA1", lambda1) ||
				!GetParam(mcParamFN, "SPECIES", "LAMBDA2", lambda2) ||
				!GetParam(mcParamFN, "SPECIES", "LAMBDA3", lambda3) ||
				!SetParam(srcParamFN, "SPECIES", "LAMBDA1", lambda1) ||
				!SetParam(srcParamFN, "SPECIES", "LAMBDA2", lambda2) ||
				!SetParam(srcParamFN, "SPECIES", "LAMBDA3", lambda3))
				return msgErrorStr(-14, "LAMBDA", srcParamFN);
			
			msgText(("MCalcs Complete, lambda1 = " + toStr(lambda1) + ", lambda2 = " + toStr(lambda2) + ", lambda3 = " + toStr(lambda3)).c_str());
		}
		else {
			//Run MCalcs generating single lambda
			result = RunMCalcs(mcParamFN);
			if (result != 0) return msgErrorStr(-15, "RunMCalcs()", toStr(result));

			//Get lamda parameter from MCalcs parameter file and set in source parameter file
			if (!GetParam(mcParamFN, "SPECIES", "LAMBDA", lambda) ||
				!SetParam(srcParamFN, "SPECIES", "LAMBDA", lambda))
				return msgErrorStr(-14, "LAMBDA", srcParamFN);
			msgText(("MCalcs Complete, lambda = " + toStr(lambda)).c_str());
		}
	}

	if (onlyMCalcs) return 0;

	////////////////////////////////////////////////////////////////////////////////////
	//Create input grids for each time step T
	////////////////////////////////////////////////////////////////////////////////////
	for (int t = 0; t < timeSteps; t++) {
		std::string tStr = toStr(t);
		std::string srHabFNt(""),
			srHHaFNt(srDir + "\\" + preFN + "SR_Habitat_Area_t" + tStr),
			srPrmFNt(srDir + "\\" + preFN + "SR_Permeability_t" + tStr),
			hrHabFNt(hrDir + "\\" + preFN + "HR_Habitat_Area_t" + tStr),
			hrPrmFNt(hrDir + "\\" + preFN + "HR_Permeability_t" + tStr),
			hrNhaFNt(hrDir + "\\" + preFN + "HR_NHA_t" + tStr),
			hrMpcFNt(hrDir + "\\" + preFN + "HR_MPC_t" + tStr),
			drHabFNt(drDir + "\\" + preFN + "DR_Habitat_Area_t" + tStr),
			drPrmFNt(drDir + "\\" + preFN + "DR_Permeability_t" + tStr),
			drNhaFNt(drDir + "\\" + preFN + "DR_NHA_t" + tStr),
			drOccFNt(drDir + "\\" + preFN + "DR_Occupancy_t" + tStr),
			drMpcFNt(drDir + "\\" + preFN + "DR_MPC_t" + tStr);

		if (!GetParam(srcParamFN, "INPUTS", "HAB_FN" + tStr, srHabFNt)) return msgErrorStr(-1, "HAB_FN" + tStr, srcParamFN);
		
		//Read in tif habitat grids. Hack job for SoS PLP project
		if (GetFileExt(srHabFNt).compare("tif") == 0) {

			if (t > 0) {
				igstream inTif;
				ogstream outFlt;

				std::cout << "Copying tif habitat grid to floating point grid\n";

				inTif.open(srHabFNt);
				if (!inTif.is_open()) return msgErrorStr(-3, srHabFNt);
				if (inTif.dataType() != 6) return msgErrorStr(-99, "Input tif habitat grid data type is not float32");

				outFlt.open(GridFN(srDir + "\\" + preFN + "SR_Habitat_t" + tStr));
				outFlt.copyHeader(inTif);
				if (!outFlt.is_open()) return msgErrorStr(-4, srDir + "\\" + preFN + "SR_Habitat_t" + tStr);

				std::unique_ptr<float[]> data = std::make_unique<float[]>(inTif.nCols());

				for (uint r = 0; r < inTif.nRows(); r++) {
					inTif.read((char*)data.get(), sizeof(float) * inTif.nCols());
					outFlt.write((const char*)data.get(), sizeof(float) * inTif.nCols());
				}

				inTif.close();
				outFlt.close();
			}
			srHabFNt = (srDir + "\\" + preFN + "SR_Habitat_t" + tStr);
		}
		if (!FileExists(GridFN(srHabFNt))) return msgErrorStr(-3, srHabFNt);
		//End of hack


		if (doSRInputs) {
			if (applyAreaToHab) {
				std::cout << "Rescaling source resolution habitat suitability grid to habitat area:\n" << srHHaFNt << "\n\n";
				if (!ApplyLambdaToGrid(srHabFNt, srHHaFNt, F2FLAMBDA(float(srCellArea) * x / srHabMax;))) return msgErrorStr(-13, srHHaFNt);
			}

			//Create perm grid (SR, T) if needed
			if (!GetParam(srcParamFN, "INPUTS", "PRM_FN" + tStr, srPrmFNt) || !FileExists(GridFN(srPrmFNt))) {
				std::cout << "Rescaling source resolution habitat suitability grid to habitat permeability:\n" << srPrmFNt << "\n\n";
				if (!CreatePrmGrid(srHabFNt, srPrmFNt, hrDistMin, hrDistMax)) return msgErrorStr(-13, srPrmFNt);
				if (!SetParam(srcParamFN, "INPUTS", "PRM_FN" + tStr, srPrmFNt)) return msgErrorStr(-14, "PRM_FN" + tStr, srcParamFN);
			}
		}
		//If source cellsize = homerange cellsize then use source habitat and permeability for homerange 
		if (doHRInputs) {
			if (hrCellSize == srCellSize) {
				std::cout << "Using source resolution habitat area and permeability grids for homerange analysis:\n" << srHHaFNt << "\n" << srPrmFNt << "\n\n";
				hrHabFNt = srHHaFNt;
				hrPrmFNt = srPrmFNt;
			}
			//Else create HR habitat and permeability and run HR context CBA
			else {
				//Check required inputs exist
				if (!FileExists(GridFN(srHabFNt)) || !FileExists(GridFN(srPrmFNt)))
					return msgErrorStr(-99, "Files required to create homerange inputs not found:\n", srHabFNt + "\n" + srPrmFNt);

				//Resample source habitat grid to homerange habitat grid, Method = Average
				std::cout << "Resampling source resolution habitat suitability grid to homerange habitat area:\n" << hrHabFNt << "\n\n";
				if (!ResampleGrid(srHabFNt, hrHabFNt, hrCellSize, 5)) return msgErrorStr(-13, hrHabFNt);
				if (!ApplyLambdaToGrid(hrHabFNt, F2FLAMBDA(float(hrCellArea) * x / srHabMax;))) return msgErrorStr(-13, hrHabFNt);

				//Resample source permeability grid to homerange permeability grid, Method = Average
				std::cout << "Resampling source resolution habitat permeability grid to homerange permeability:\n" << hrPrmFNt << "\n\n";
				if (!ResampleGrid(srPrmFNt, tmpFN, hrHabFNt, 5)) return msgErrorStr(-13, tmpFN);
				if (!CreatePrmGrid(tmpFN, hrPrmFNt, hrDistMin, hrDistMax)) return msgErrorStr(-13, hrPrmFNt);
				CopyFS(tmpFN, hrPrmFNt + "_Temp");
				DeleteGrid(tmpFN);
			}
		}

		//Perform context CBA using homerange habitat and permeability grids to create homerange NHA grid
		if (doHomerange) {
			if (doHRIterations) {
				//Write filenames to homerange REMP CBA parameter file
				if (!SetParam(hrParamFN, "INPUTS", "HAB_FN0", hrHabFNt) ||
					!SetParam(hrParamFN, "INPUTS", "PRM_FN0", hrPrmFNt) ||
					!SetParam(hrParamFN, "OUTPUTS", "MPC_FN0", hrMpcFNt) ||
					!SetParam(hrParamFN, "WINDOW", "LOOKUP_FN", hrLookupFN) ||
					!SetParam(hrParamFN, "WINDOW", "TRANSLATE_FN", hrTranslateFN) ||
					!SetParam(hrParamFN, "SPECIES", "LAMBDA", 1.0) || //REMP_TS_CBA() wants LAMBDA but it's not needed here, make optional
					!SetParam(hrParamFN, "SPECIES", "#NOTE", "LAMBDA is not used here but parameter needed by REMP CBA") || 
					!SetParam(hrParamFN, "ANALYSIS", "ITERATIONS", hrIterations) ||
					!SetParam(hrParamFN, "ANALYSIS", "MAXPETALPERM", useMaxPetalPerm) ||
					!SetParam(hrParamFN, "ANALYSIS", "DEBUG", rempDebug))
					return msgErrorStr(-14, "parameters", hrParamFN);
				
				std::cout << "Performing homerange MPC analysis using parameters from:\n" << hrParamFN << "\n\n";
				result = OccupancyCBA_TS(hrParamFN.c_str());
				if (result != 0) return msgErrorStr(-15, "REMP_TS_CBA()", toStr(result));
				hrNhaFNt = hrMpcFNt;
			}
			else {
				//Check required inputs exist
				if (!FileExists(GridFN(hrHabFNt)) || !FileExists(GridFN(hrPrmFNt)))
					return msgErrorStr(-99, "Files required to process homerange not found:\n", hrHabFNt + "\n" + hrPrmFNt);

				if (!SetParam(hrParamFN, "INPUTS", "HAB_FN", hrHabFNt) ||
					!SetParam(hrParamFN, "INPUTS", "PRM_FN", hrPrmFNt) ||
					!SetParam(hrParamFN, "OUTPUTS", "NHA_FN", hrNhaFNt) ||
					!SetParam(hrParamFN, "WINDOW", "LOOKUP_FN", hrLookupFN) ||
					!SetParam(hrParamFN, "WINDOW", "TRANSLATE_FN", hrTranslateFN))
					return msgErrorStr(-14, "parameters", hrParamFN);

				std::cout << "Performing homerange NHA analysis using parameters from:\n" << hrParamFN << "\n\n";
				result = ContextCBA(hrParamFN.c_str());
				if (result != 0) return msgErrorStr(-15, "ContextCBA()", toStr(result));
			}
		}

		if (doDRInputs) {
			//Check required inputs exist
			if (!FileExists(GridFN(srHabFNt)) || !FileExists(GridFN(srPrmFNt)))
				return msgErrorStr(-99, "Files required to create dispersal inputs not found:\n", srHabFNt + "\n" + srPrmFNt);
			if (!GetParam(srcParamFN, "INPUTS", "NHA_FN" + tStr, hrNhaFNt) && !FileExists(GridFN(hrNhaFNt)))
				return msgErrorStr(-99, "Files required to create dispersal inputs not found:\n", hrNhaFNt);

			//Resample source habitat grid to dispersal habitat grid, Method = Average
			std::cout << "Resampling source resolution habitat suitability grid to dispersal habitat area:\n" << drHabFNt << "\n\n";
			if (!ResampleGrid(srHabFNt, drHabFNt, drCellSize, 5)) return msgErrorStr(-13, drHabFNt);
			if (!ApplyLambdaToGrid(drHabFNt, F2FLAMBDA(float(drCellArea) * x / srHabMax;))) return msgErrorStr(-13, drHabFNt);

			//Resample source permeability grid to dispersal permeability grid, Method = Average
			std::cout << "Resampling source resolution habitat permeability grid to dispersal permeability:\n" << drPrmFNt << "\n\n";
			if (!ResampleGrid(srPrmFNt, tmpFN, drHabFNt, 5)) return msgErrorStr(-13, tmpFN);
			if (!CreatePrmGrid(tmpFN, drPrmFNt, drDistMin, drDistMax)) return msgErrorStr(-13, hrPrmFNt);
			CopyFS(tmpFN, drPrmFNt + "_Temp");
			DeleteGrid(tmpFN);

			//Resample homerange (or source) NHA grid to dispersal NHA grid, Method = Average
			std::cout << "Resampling homerange NHA grid to dispersal NHA:\n" << drNhaFNt << "\n\n";
			if (!ResampleGrid(hrNhaFNt, drNhaFNt, drHabFNt, 5)) return msgErrorStr(-13, drNhaFNt);

			//Write filenames to dispersal REMP CBA parameter file
			if (!SetParam(drParamFN, "INPUTS", "HAB_FN" + tStr, drHabFNt) ||
				!SetParam(drParamFN, "INPUTS", "PRM_FN" + tStr, drPrmFNt) ||
				!SetParam(drParamFN, "INPUTS", "NHA_FN" + tStr, drNhaFNt) ||
				!SetParam(drParamFN, "OUTPUTS", "OCC_FN" + tStr, drOccFNt) ||
				!SetParam(drParamFN, "OUTPUTS", "MPC_FN" + tStr, drMpcFNt))
				return msgErrorStr(-14, "parameters", drParamFN);
		}
	}

	//Perform remp timeseries CBA using habitat grid(DR, T), permeability grid(DR, T), NHA grid(DR, T) to create Occ grid(DR, T) and MPC grid(DR, T)
	if (doDispersal) {	
		//Check required inputs exist
		for (int t = 0; t < timeSteps; t++) {
			std::string tStr = toStr(t),
				drHabFNt(drDir + "\\" + preFN + "DR_Habitat_Area_t" + tStr),
				drPrmFNt(drDir + "\\" + preFN + "DR_Permeability_t" + tStr),
				drNhaFNt(drDir + "\\" + preFN + "DR_NHA_t" + tStr);
			if (!FileExists(GridFN(drHabFNt)) || !FileExists(GridFN(drPrmFNt)) || !FileExists(GridFN(drNhaFNt)))
				return msgErrorStr(-99, "Files required to process dispersal not found:\n", drHabFNt + "\n" + drPrmFNt + "\n" + drNhaFNt);
		}

		if (!SetParam(drParamFN, "WINDOW", "LOOKUP_FN", drLookupFN) ||
			!SetParam(drParamFN, "WINDOW", "TRANSLATE_FN", drTranslateFN) ||
			!SetParam(drParamFN, "SPECIES", "LAMBDA", doTriLambda ? lambda2 : lambda) ||
			!SetParam(drParamFN, "ANALYSIS", "TIMESTEPS", timeSteps) ||
			!SetParam(drParamFN, "ANALYSIS", "ITERATIONS", iterations) ||
			!SetParam(drParamFN, "ANALYSIS", "MAXPETALPERM", useMaxPetalPerm) ||
			!SetParam(drParamFN, "ANALYSIS", "COUPLED", coupled) ||
			!SetParam(drParamFN, "ANALYSIS", "DEBUG", rempDebug))
			return msgErrorStr(-14, "parameters", drParamFN);
		
		std::cout << "Performing dispersal analysis using parameters from:\n" << drParamFN << "\n\n";
		result = OccupancyCBA_TS(drParamFN.c_str());
		if (result != 0) return msgErrorStr(-15, "REMP_TS_CBA()", toStr(result));

		//	Do Pu (DR_Occupancy_t) = Si / (Si + lambda3 / NHA)  ==  Si * NHA / (Si * NHA + lambda3) 
		if (doTriLambda) {
			for (int t = 0; t < timeSteps; t++) {
				std::string tStr = toStr(t);
				std::string SiFNt(drDir + "\\" + preFN + "DR_Occupancy_t" + tStr + "_Si");
				std::string nhaFNt(drDir + "\\" + preFN + "DR_NHA_t" + tStr);
				std::string drOccFNt(drDir + "\\" + preFN + "DR_Occupancy_t" + tStr);

				if (!FileExists(GridFN(SiFNt)) || !FileExists(GridFN(nhaFNt)))
					return msgErrorStr(-99, "Files required to create Pu grid not found:\n", SiFNt + "\n" + nhaFNt);
				std::cout << "Creating dispersal resolution occupancy (Pu) grid:\n" << drOccFNt << "\n\n";
				if (!ApplyLambdaToTwoGrids(SiFNt, nhaFNt, drOccFNt, FF2FLAMBDA(float((double(a) * double(b)) / (double(a) * double(b) + lambda3));))) return msgErrorStr(-13, drOccFNt);
			}
		}
	}

	//For each time step T, resample Pu grid to SR/HR then apply 1 - (1 - Pu grid)^habitat grid to create Pi grid
	if (doFinalPi) {
		float habDen;
		std::string habFNt("");
		for (int t = 0; t < timeSteps; t++) {
			std::string tStr = toStr(t);
			habDen = srHabMax;
			if (doFinalStepAtHR && hrCellSize != srCellSize) {
				habFNt = hrDir + "\\" + preFN + "HR_Habitat_Area_t" + tStr;
				habDen = float(hrCellArea);
			}
			else if (applyAreaToHab) {
				habFNt = srDir + "\\" + preFN + "SR_Habitat_Area_t" + tStr;
				habDen = float(srCellArea);
			}
			else if (!GetParam(srcParamFN, "INPUTS", "HAB_FN" + tStr, habFNt)) 
				return msgErrorStr(-1, "HAB_FN" + tStr, srcParamFN);

			if (GetFileExt(habFNt).compare("tif") == 0) {
				habFNt = (srDir + "\\" + preFN + "SR_Habitat_t" + tStr);
			}
			
			std::string drOccFNt(drDir + "\\" + preFN + "DR_Occupancy_t" + tStr);
			
			//Following changes to write final grids as tifs until full gstream implementation 
			//std::string PuFNt(ocDir + "\\" + preFN + "Pu_t" + tStr);
			//std::string PiFNt(ocDir + "\\" + preFN + "Pi_t" + tStr);
			std::string PuFNt(ocDir + "\\" + preFN + "Pu_t" + tStr + ".tif");
			std::string PiFNt(ocDir + "\\" + preFN + "Pi_t" + tStr + ".tif");

			if (!FileExists(GridFN(habFNt)) || !FileExists(GridFN(drOccFNt)))
				return msgErrorStr(-99, "Files required to create Pi grid not found:\n", habFNt + "\n" + drOccFNt);
			
			std::cout << "Resampling dispersal resolution occupancy (Pu) grid to source resolution:\n" << PuFNt << "\n\n";
			if (!ResampleGrid(drOccFNt, PuFNt, habFNt, 1)) return msgErrorStr(-13, PuFNt);

			std::cout << "Creating source resolution Pi grid:\n" << PiFNt << "\n\n";
			if (!ApplyLambdaToTwoGStreams(PuFNt, habFNt, PiFNt, FF2FLAMBDA(1.0f - std::pow(1.0f - a, b / habDen);))) return msgErrorStr(-13, PiFNt);

		}
	}

	//Cleanup
	std::cout << "performing any clean up of folders\n\n";
	if (cleanAll || cleanMC) fileSys::remove_all(mcDir);
	if (cleanAll || cleanSR) fileSys::remove_all(srDir);
	if (cleanAll || cleanHR) fileSys::remove_all(hrDir);
	if (cleanAll || cleanDR) fileSys::remove_all(drDir);

	//Complete!
	msgText("RunREMP() Complete!");
	return result;
}