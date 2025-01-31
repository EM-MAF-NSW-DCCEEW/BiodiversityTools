/*
Profiler.h - Profiler class decleration
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
#ifndef PROFILER_H
#define PROFILER_H

#ifdef _WIN32
#	include <windows.h>
#else	//LINUX
#	include <chrono>
#endif


class Profiler
{
private:
#ifdef _WIN32
	LARGE_INTEGER freq;
	LARGE_INTEGER start;
	LARGE_INTEGER elapsed;
#else	//LINUX
	std::chrono::time_point<std::chrono::steady_clock> start;
	std::chrono::time_point<std::chrono::steady_clock> elapsed;
#endif
	long long totalTime;
	long long elapsedTime;
	long long iterations;
	long long metric;
	bool isRunning;

public:
	Profiler();
	Profiler(int t);
	~Profiler(void);
	void Start();
	long long Elapsed();
	long long Stop();
	long long Total();
	long long Count();
	double Avg();
};
//#else use GNU counters if building for *nix
#endif