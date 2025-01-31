/*
SpatialLinksLib.h - entrypoint functions for Spatial Links library
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

//Notes
//TODO

#pragma once
#ifndef SPATIAL_LINKS_LIB_H
#define SPATIAL_LINKS_LIB_H

#include <string>

#ifdef _WIN32
#   ifdef LINKS_EXPORTS
#       define LINKS_API __declspec(dllexport)
#   else
#       define LINKS_API __declspec(dllimport)
#   endif
#   define STDCALL __stdcall
#else //LINUX
#   define LINKS_API 
#   define STDCALL 
#endif

//Message function pointer types
typedef void(STDCALL *msgTextFP)(const char *msg);
typedef void(STDCALL *msgProgressFP)(const char *msg, int val);

//Entrypoint functions for Spatial Links library
LINKS_API int LinksCompleteSampling(char *paramFN);
LINKS_API int LinksCompleteSampling(char *paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB);

LINKS_API int LinksCompleteSamplingAltDstHab(char *paramFN);
LINKS_API int LinksCompleteSamplingAltDstHab(char *paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB);

LINKS_API int LinksRandomSampling(char *paramFN);
LINKS_API int LinksRandomSampling(char *paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB);

#endif
