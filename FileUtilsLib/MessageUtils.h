/*
MessageUtils.h - Message utility function declerations
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
#ifndef MESSAGE_UTILS_H
#define MESSAGE_UTILS_H
#include <iostream>
#include <string>

#ifdef _WIN32
#   define STDCALL __stdcall
#else //LINUX
#   define STDCALL 
#endif

//case 0:	msg = "No error"; break;
//case -1:	msg = "Error: Parameter " + arg1 + " not found in " + arg2; break;
//case -2:	msg = "Error: Filename not valid: " + arg1; break;
//case -3:	msg = "Error: File not readable: " + arg1; break;
//case -4:	msg = "Error: File not writable: " + arg1; break;
//case -5:	msg = "Error: Header file not readable: " + arg1; break;
//case -6:	msg = "Error: Header file not writable: " + arg1; break;
//case -7:	msg = "Error: Header files don't match: " + arg1 + " " + arg2; break;
//case -8:	msg = "Error: Unable to copy file " + arg1 + " to " + arg2; break;
//case -9:	msg = "Error: Unable to create petals " + arg1; break;
//case -10:	msg = "Error: Unable to access CUDA device"; break;
//case -11:	msg = "Error: Exception caught: " + arg1; break;
//case -12:	msg = "Error: Unable to allocate filenames"; break;
//case -13:	msg = "Error: Unable to create file: " + arg1; break;
//case -14:	msg = "Error: Unable to set " + arg1 + " in " + arg2; break;
//case -15:	msg = "Error: Function " + arg1 + " returned non-zero error code: " + arg2; break;
int msgError(int errorCode, const char *arg1 = "", const char *arg2 = "", const char *arg3 = "");
//case 0:	msg = "No error"; break;
//case -1:	msg = "Error: Parameter " + arg1 + " not found in " + arg2; break;
//case -2:	msg = "Error: Filename not valid: " + arg1; break;
//case -3:	msg = "Error: File not readable: " + arg1; break;
//case -4:	msg = "Error: File not writable: " + arg1; break;
//case -5:	msg = "Error: Header file not readable: " + arg1; break;
//case -6:	msg = "Error: Header file not writable: " + arg1; break;
//case -7:	msg = "Error: Header files don't match: " + arg1 + " " + arg2; break;
//case -8:	msg = "Error: Unable to copy file " + arg1 + " to " + arg2; break;
//case -9:	msg = "Error: Unable to create petals " + arg1; break;
//case -10:	msg = "Error: Unable to access CUDA device"; break;
//case -11:	msg = "Error: Exception caught: " + arg1; break;
//case -12:	msg = "Error: Unable to allocate filenames"; break;
//case -13:	msg = "Error: Unable to create file: " + arg1; break;
//case -14:	msg = "Error: Unable to set " + arg1 + " in " + arg2; break;
//case -15:	msg = "Error: Function " + arg1 + " returned non-zero error code: " + arg2; break;
int msgErrorStr(int errorCode, const std::string & arg1 = "", const std::string & arg2 = "", const std::string & arg3 = "");

template<typename T> inline std::string toStr(T v) { return std::to_string(v); }

//text, update and progress message function pointer types 
typedef void(STDCALL *msgTextFP)(const char *msg);
//typedef void(STDCALL *msgTextFP)(const std::string &msg);
typedef void(STDCALL *msgProgressFP)(const char *msg, int val);
typedef int (STDCALL *msgGDALProgressFP)(double dfComplete, const char *pszMessage, void *pProgressArg);

typedef void(STDCALL *msgStringFP)(const std::string &str);

//Declare extern message function pointers
extern msgTextFP msgText;
extern msgProgressFP msgProgress;
extern msgGDALProgressFP msgGDALProgress;

extern msgStringFP msgString;

#endif