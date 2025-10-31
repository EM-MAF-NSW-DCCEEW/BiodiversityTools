/*
SpatialContextLib.h - Exported functions for SpatialContextLib
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

//NOTE Do NOT compile with --use_fast_math in CUDA 9.2 for c5.2_sm5.2 maybe others, causes non-reproducable errors

#pragma once
#ifndef _SPATIAL_CONTEXT_LIB_H_
#define _SPATIAL_CONTEXT_LIB_H_
#include <string>
#include <vector>


#ifdef _WIN32
#   ifdef CUDA_CBA_EXPORTS
#       define CUDA_CBA_API __declspec(dllexport)
#   else
#       define CUDA_CBA_API __declspec(dllimport)
#   endif
#   define STDCALL __stdcall
#else //LINUX
#   define CUDA_CBA_API 
#   define STDCALL 
#endif

//Message function pointer types
typedef void(STDCALL *msgTextFP)(const char *msg);
typedef void(STDCALL *msgProgressFP)(const char *msg, int val);

//Exported functions
//Context CBA
CUDA_CBA_API int ContextCBA(const char *paramFN);
CUDA_CBA_API int ContextCBA(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB);

//Benfits CBA
CUDA_CBA_API int BenefitsCBA(const char *paramFN);
CUDA_CBA_API int BenefitsCBA(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB);

//EHA CBA
CUDA_CBA_API int EHA_CBA(const char *paramFN);
CUDA_CBA_API int EHA_CBA(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB);

//REMP CBA
CUDA_CBA_API int OccupancyCBA_TS(const char *paramFN);
CUDA_CBA_API int OccupancyCBA_TS(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB);

////Resilience CBA
//CUDA_CBA_API int ResilienceCBA(const char *paramFN);
//CUDA_CBA_API int ResilienceCBA(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB);
//
////Resilience CBA 2
//CUDA_CBA_API int ResilienceCBA2(const char *paramFN);
//CUDA_CBA_API int ResilienceCBA2(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB);
//
////Resilience Denom
//CUDA_CBA_API int ResilienceDenom(const char *paramFN);
//CUDA_CBA_API int ResilienceDenom(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB);
//
////GDM Similarity to Class
//CUDA_CBA_API int GDMSimToClass(const char *paramFN);
//CUDA_CBA_API int GDMSimToClass(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB);
//
////GDM Similarity to Pixel
//CUDA_CBA_API int GDMSimToPixel(const char *paramFN);
//CUDA_CBA_API int GDMSimToPixel(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB);
//
////GDM Similarity to Species Observations
//CUDA_CBA_API int GDMSimToSpcObs(const char *paramFN);
//CUDA_CBA_API int GDMSimToSpcObs(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB);

//
////Uniqueness
CUDA_CBA_API int Uniqueness(const char* paramFN);
CUDA_CBA_API int Uniqueness(const char* paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB);

//
////Continuous MBV
CUDA_CBA_API int ContinuousMBV(const char *paramFN);
CUDA_CBA_API int ContinuousMBV(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB);

//
////Continuous BDI
//CUDA_CBA_API int ContinuousBDI(const char *paramFN);
//CUDA_CBA_API int ContinuousBDI(const char * paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB);
//
//CUDA_CBA_API int IterateSumNHab(const char *paramFN);

//Create Petals and Segments
CUDA_CBA_API int CreatePetals(const char *outFN, int srcSize, int dstSize, float zoneRatio, bool createAllFiles);

CUDA_CBA_API int GetPetalFiles(const std::string & paramFN, const std::string & paramSection, const std::string workDir, std::string & lookupFN, std::string & translateFN);
CUDA_CBA_API int GetPetalFiles(const std::string & paramFN, const std::string & paramSection, std::string & lookupFN, std::string & translateFN);

CUDA_CBA_API int CreateSegments(int radCells, int nRings, int nSlices, std::vector<double> radList, const char *outFN);
CUDA_CBA_API int CreateSegments(int radCells, int nRings, int nSlices, double zoneRatio, const char *outFN);

//Search window functions
CUDA_CBA_API int CreateSearchWindow(const std::string& paramFN);
CUDA_CBA_API int CreateSearchWindow(const std::string& paramFN, const std::string& paramSection);


#endif

