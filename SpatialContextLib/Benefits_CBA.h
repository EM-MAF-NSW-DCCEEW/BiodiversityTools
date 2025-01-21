/*
Benefits_CBA.h - CBA functions for performing BFT benefits analysis
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

#pragma once
#ifndef _BENEFITS_CBA_H_
#define _BENEFITS_CBA_H_
#include "Parameters.h"

//Benefits CBA Functions
//int BenefitsCBA(char *pFN);
int BenefitsCBA(const char *paramFN, CBAParams &params);
int BenefitsCBA(const std::string & paramFN, CBAParams &params);

#endif