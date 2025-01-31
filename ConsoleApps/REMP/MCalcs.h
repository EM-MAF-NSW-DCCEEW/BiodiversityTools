/*
MCalcs.h - TODO Description
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

#pragma once
#ifndef MCALCS_H
#define MCALCS_H

//Creates MCalcs input grids then runs REMP CBA to calculate lambda and c
int RunMCalcs(const char *paramFN);
int RunMCalcs(const std::string & paramFN);
int RunMCalcsTriLambda(const std::string & paramFN);
#endif