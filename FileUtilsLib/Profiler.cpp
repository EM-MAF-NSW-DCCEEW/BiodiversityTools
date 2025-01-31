/*
Profiler.cpp - Profiler class member definitions
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


#include "Profiler.h"
#include <iostream>

Profiler::Profiler()
{
	//Default microseconds
	Profiler(1000000);
}

Profiler::Profiler(int t)
{
#ifdef _WIN32
	QueryPerformanceFrequency(&freq);
#else	//LINUX
	//std::cout << "Linux Profiler not implemented" << std::endl;
#endif

	metric = (long long)t;
	totalTime = 0LL;
	elapsedTime = 0LL;
	iterations = 0LL;
	isRunning = false;
}

Profiler::~Profiler(void){};

void Profiler::Start()
{
	if (isRunning) Stop();
	isRunning = true;
	iterations++;
#ifdef _WIN32
	QueryPerformanceCounter(&start);
#else	//LINUX
	//std::cout << "Linux Profiler not implemented" << std::endl;
	start = std::chrono::steady_clock::now();
#endif
}

long long Profiler::Elapsed()
{
#ifdef _WIN32
	QueryPerformanceCounter(&elapsed);
	elapsedTime = (elapsed.QuadPart - start.QuadPart) * metric / freq.QuadPart;
#else	//LINUX
	//std::cout << "Linux Profiler not implemented" << std::endl;
	elapsed = std::chrono::steady_clock::now();
	elapsedTime = std::chrono::duration_cast<std::chrono::microseconds>(elapsed - start).count();
#endif
	return elapsedTime;
}

long long Profiler::Stop()
{
#ifdef _WIN32
	QueryPerformanceCounter(&elapsed);
	elapsedTime = (elapsed.QuadPart - start.QuadPart) * metric / freq.QuadPart;
#else	//LINUX
	//std::cout << "Linux Profiler not implemented" << std::endl;
	elapsed = std::chrono::steady_clock::now();
	elapsedTime = std::chrono::duration_cast<std::chrono::microseconds>(elapsed - start).count(); 
#endif
	totalTime += elapsedTime;
	isRunning = false;
	return elapsedTime;
}

long long Profiler::Total(){ return totalTime; }
long long  Profiler::Count() { return iterations; }
double Profiler::Avg(){ return (iterations == 0 ? 0 : (double)totalTime / (double)iterations); }
