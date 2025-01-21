/*
REMP_CBA.h - CBA functions for performing REMP analysis
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
#ifndef _OCCUPANCY_CBA_H_
#define _OCCUPANCY_CBA_H_
#include "Parameters.h"

//REMP CBA Functions
int OccupancyCBA_TS(const char *paramFN, CBAParams &params);
int OccupancyCBA_TS(const std::string & paramFN, CBAParams &params);
#endif
