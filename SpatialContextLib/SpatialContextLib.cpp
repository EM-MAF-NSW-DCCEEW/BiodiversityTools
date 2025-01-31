/*
SpatialContextLib.cpp - Exported functions for SpatialContextLib
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

//NOTE Do NOT use --use_fast_math in CUDA 9.2 for c5.2_sm5.2 causes non-reproducable errors

//Include all CBA function header files
#include "Petals.h"
#include "FileUtilsLib.h"
#include "SpatialContextLib.h"
#include "Context_CBA.h"
#include "Benefits_CBA.h"
#include "EHA_CBA.h"
#include "Occupancy_CBA.h"

//#include "Resilience_CBA.h"
//#include "Resilience_CBA2.h"
//#include "ResilienceDenominator.h"
//#include "GDMSimToClass.h"
//#include "GDMSimToPixel.h"
//#include "GDMSimToSpcObs.h"
//#include "ContinuousMBV.h"
//#include "ContinuousBDI.h"
//#include "IterateSumNHab.h"

//Exported functions

//Context CBA exported functions
CUDA_CBA_API int ContextCBA(const char *paramFN) {
	CBAParams params;
	return ContextCBA(paramFN, params);
}

CUDA_CBA_API int ContextCBA(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB) {
	CBAParams params;
	msgText = msgTextCB;
	msgProgress = msgUpdateCB;
	return ContextCBA(paramFN, params);
}

//Benefits CBA exported functions
CUDA_CBA_API int BenefitsCBA(const char *paramFN) {
	CBAParams params;
	return BenefitsCBA(paramFN, params);
}

CUDA_CBA_API int BenefitsCBA(const char *paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB) {
	CBAParams params;
	msgText = msgTextCB;
	msgProgress = msgUpdateCB;
	return BenefitsCBA(paramFN, params);
}

//EHA CBA exported functions
CUDA_CBA_API int EHA_CBA(const char *paramFN) {
	CBAParams params;
	return EHA_CBA(paramFN, params);
}

CUDA_CBA_API int EHA_CBA(const char *paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB) {
	CBAParams params;
	msgText = msgTextCB;
	msgProgress = msgUpdateCB;
	return EHA_CBA(paramFN, params);
}

//REMP CBA exported functions
CUDA_CBA_API int OccupancyCBA_TS(const char *paramFN) {
	CBAParams params;
	return OccupancyCBA_TS(paramFN, params);
}

CUDA_CBA_API int OccupancyCBA_TS(const char *paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB) {
	CBAParams params;
	msgText = msgTextCB;
	msgProgress = msgUpdateCB;
	return OccupancyCBA_TS(paramFN, params);
}

////Resilience CBA exported functions
//CUDA_CBA_API int ResilienceCBA(const char *paramFN) {
//	CBAParams params;
//	return ResilienceCBA(paramFN, params);
//}
//
//CUDA_CBA_API int ResilienceCBA(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB) {
//	CBAParams params;
//	msgText = msgTextCB;
//	msgProgress = msgUpdateCB;
//	return ResilienceCBA(paramFN, params);
//}
//
////Resilience CBA2 exported functions
//CUDA_CBA_API int ResilienceCBA2(const char *paramFN) {
//	CBAParams params;
//	return ResilienceCBA2(paramFN, params);
//}
//
//CUDA_CBA_API int ResilienceCBA2(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB) {
//	CBAParams params;
//	msgText = msgTextCB;
//	msgProgress = msgUpdateCB;
//	return ResilienceCBA2(paramFN, params);
//}
//
////Resilience Denominator exported functions
//CUDA_CBA_API int ResilienceDenom(const char *paramFN) {
//	CBAParams params;
//	return ResilienceDenom(paramFN, params);
//}
//
//CUDA_CBA_API int ResilienceDenom(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB) {
//	CBAParams params;
//	msgText = msgTextCB;
//	msgProgress = msgUpdateCB;
//	return ResilienceDenom(paramFN, params);
//}
//
////GDM Similarity to class exported functions
//CUDA_CBA_API int GDMSimToClass(const char *paramFN) {
//	CBAParams params;
//	return GDMSimToClass(paramFN, params);
//}
//
//CUDA_CBA_API int GDMSimToClass(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB) {
//	CBAParams params;
//	msgText = msgTextCB;
//	msgProgress = msgUpdateCB;
//	return GDMSimToClass(paramFN, params);
//}
//
//// GDM Similarity to pixel exported functions
//CUDA_CBA_API int GDMSimToPixel(const char *paramFN) {
//	CBAParams params;
//	return GDMSimToPixel(paramFN, params);
//}
//
//CUDA_CBA_API int GDMSimToPixel(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB) {
//	CBAParams params;
//	msgText = msgTextCB;
//	msgProgress = msgUpdateCB;
//	return GDMSimToPixel(paramFN, params);
//}
//
////GDM Similarity to species observations exported functions
//CUDA_CBA_API int GDMSimToSpcObs(const char *paramFN) {
//	CBAParams params;
//	return GDMSimToSpcObs(paramFN, params);
//}
//
//CUDA_CBA_API int GDMSimToSpcObs(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB) {
//	CBAParams params;
//	msgText = msgTextCB;
//	msgProgress = msgUpdateCB;
//	return GDMSimToSpcObs(paramFN, params);
//}
//
////Continuous MBV exported functions
//CUDA_CBA_API int ContinuousMBV(const char *paramFN) {
//	CBAParams params;
//	return ContinuousMBV(paramFN, params);
//}
//
//CUDA_CBA_API int ContinuousMBV(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB) {
//	CBAParams params;
//	msgText = msgTextCB;
//	msgProgress = msgUpdateCB;
//	return ContinuousMBV(paramFN, params);
//}
//
////Continuous BDI exported functions
//CUDA_CBA_API int ContinuousBDI(const char *paramFN) {
//	CBAParams params;
//	return ContinuousBDI(paramFN, params);
//}
//
//CUDA_CBA_API int ContinuousBDI(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB) {
//	CBAParams params;
//	msgText = msgTextCB;
//	msgProgress = msgUpdateCB;
//	return ContinuousBDI(paramFN, params);
//}
//
////IterateSumNHab exported functions
//CUDA_CBA_API int IterateSumNHab(const char *paramFN) {
//	CBAParams params;
//	return IterateSumNHab(paramFN, params);
//}

//Create Petals and Segments exported functions
CUDA_CBA_API int CreatePetals(const char *outFN, int srcSize, int dstSize, float zoneRatio, bool createAllFiles) {
	CBAParams params;
	return CreatePetals(outFN, srcSize, dstSize, zoneRatio, createAllFiles, params);
}

CUDA_CBA_API int GetPetalFiles(const std::string & paramFN, const std::string & paramSection, const std::string workDir, std::string & lookupFN, std::string & translateFN) {
	CBAParams params;
	return GetPetalFiles(paramFN, paramSection, workDir, lookupFN, translateFN, params);
}

CUDA_CBA_API int GetPetalFiles(const std::string & paramFN, const std::string & paramSection, std::string & lookupFN, std::string & translateFN) {
	CBAParams params;
	std::string workDir;
	GetFilesFolder(paramFN, workDir);
	return GetPetalFiles(paramFN, paramSection, workDir, lookupFN, translateFN, params);
}

CUDA_CBA_API int CreateSegments(int radCells, int nRings, int nSlices, std::vector<double> radList, const char *outFN) {
	CBAParams params;
	return CreateSegments(radCells, nRings, nSlices, radList, outFN, params);
}

CUDA_CBA_API int CreateSegments(int radCells, int nRings, int nSlices, double zoneRatio, const char *outFN) {
	CBAParams params;
	return CreateSegments(radCells, nRings, nSlices, zoneRatio, std::string(outFN), params);
}

//Exported External C functions
extern "C" {
	//Context CBA externalC functions
	CUDA_CBA_API int ContextCBA_ExtC(const char *paramFN) {
		return ContextCBA(paramFN);
	}

	CUDA_CBA_API int ContextCBA_ExtC_CB(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB) {
		return ContextCBA(paramFN, msgTextCB, msgUpdateCB);
	}

	//Benefits CBA externalC functions
	CUDA_CBA_API int BenefitsCBA_ExtC(const char *paramFN) {
		return BenefitsCBA(paramFN);
	}

	CUDA_CBA_API int BenefitsCBA_ExtC_CB(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB) {
		return BenefitsCBA(paramFN, msgTextCB, msgUpdateCB);
	}

	//EHA CBA externalC functions
	CUDA_CBA_API int EHA_CBA_ExtC(const char *paramFN) {
		return EHA_CBA(paramFN);
	}

	CUDA_CBA_API int EHA_CBA_ExtC_CB(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB) {
		return EHA_CBA(paramFN, msgTextCB, msgUpdateCB);
	}

	//REMP CBA externalC functions
	CUDA_CBA_API int REMP_TS_CBA_ExtC(const char *paramFN) {
		return OccupancyCBA_TS(paramFN);
	}

	CUDA_CBA_API int REMP_TS_CBA_ExtC_CB(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB) {
		return OccupancyCBA_TS(paramFN, msgTextCB, msgUpdateCB);
	}

	////Resilience CBA externalC functions
	//CUDA_CBA_API int ResilienceCBA_ExtC(const char *paramFN) {
	//	return ResilienceCBA(paramFN);
	//}

	//CUDA_CBA_API int ResilienceCBA_ExtC_CB(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB) {
	//	return ResilienceCBA(paramFN, msgTextCB, msgUpdateCB);
	//}

	//Create Petals externalC functions
	CUDA_CBA_API int CreatePetals_ExtC(const char *outFN, int srcSize, int dstSize, float zoneRatio, bool createAllFiles) {
		return CreatePetals(outFN, srcSize, dstSize, zoneRatio, createAllFiles);
	}
}