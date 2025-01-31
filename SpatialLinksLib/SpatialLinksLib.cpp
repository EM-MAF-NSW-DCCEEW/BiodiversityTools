/*
SpatialLinksLib.cpp - entrypoint functions for Spatial Links library 
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
//Standardise parameter file format

//Include all Links function header files
#include "SpatialLinksLib.h"
#include "LinksParams.h"
#include "LinksCompleteSampling.h"
#include "LinksCompleteSamplingAltDstHab.h"
#include "LinksRandomSampling.h"

//Entrypoint functions for Spatial Links library
LINKS_API int LinksCompleteSampling(char *paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB) {
	msgText = msgTextCB;
	msgProgress = msgUpdateCB;
	return LinksCompleteSampling(paramFN);
}

LINKS_API int LinksCompleteSampling(char *paramFN) {
	LinksParams linksParams;
	msgText("Running Links with complete sampling");
	LinksCompleteSampling(paramFN, linksParams);
	return 0;
}

LINKS_API int LinksCompleteSamplingAltDstHab(char *paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB) {
	msgText = msgTextCB;
	msgProgress = msgUpdateCB;
	return LinksCompleteSamplingAltDstHab(paramFN);
}

LINKS_API int LinksCompleteSamplingAltDstHab(char *paramFN) {
	LinksParams linksParams;
	msgText("Running Links with complete sampling and alternate destination habitat");
	LinksCompleteSamplingAltDstHab(paramFN, linksParams);
	return 0;
}

LINKS_API int LinksRandomSampling(char *paramFN, msgTextFP msgTextCB, msgProgressFP msgUpdateCB) {
	msgText = msgTextCB;
	msgProgress = msgUpdateCB;
	return LinksRandomSampling(paramFN);
}

LINKS_API int LinksRandomSampling(char *paramFN) {
	LinksParams linksParams;
	msgText("Running Links with random sampling");
	LinksRandomSampling(paramFN, linksParams);
	return 0;
}
