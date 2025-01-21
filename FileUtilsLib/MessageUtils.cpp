/*
MessageUtils.cpp - Message utility function declerations
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


#include <algorithm>
#include "MessageUtils.h"
#pragma warning(disable : 4996)

int msgError(int errorCode, const char *arg1, const char *arg2, const char *arg3)
{
	return msgErrorStr(errorCode, std::string(arg1), std::string(arg2), std::string(arg3));
}
int msgErrorStr(int errorCode, const std::string & arg1, const std::string & arg2, const std::string & arg3)
{
	std::string msg("");
	switch (errorCode) {
	case 0:		msg = "No error"; break;
	case -1:	msg = "Error: Unable to get " + arg1 + " from " + arg2; break;
	case -2:	msg = "Error: Filename not valid: " + arg1; break;
	case -3:	msg = "Error: File not readable: " + arg1; break;
	case -4:	msg = "Error: File not writable: " + arg1; break;
	case -5:	msg = "Error: Header file not readable: " + arg1; break;
	case -6:	msg = "Error: Header file not writable: " + arg1; break;
	case -7:	msg = "Error: Header files don't match: " + arg1 + " " + arg2; break;
	case -8:	msg = "Error: Unable to copy file " + arg1 + " to " + arg2; break;
	case -9:	msg = "Error: Unable to create petals " + arg1; break;
	case -10:	msg = "Error: Unable to access CUDA device"; break;
	case -11:	msg = "Error: Exception caught: " + arg1; break;
	case -12:	msg = "Error: Unable to allocate filenames"; break;
	case -13:	msg = "Error: Unable to create file: " + arg1; break;
	case -14:	msg = "Error: Unable to set " + arg1 + " in " + arg2; break;
	case -15:	msg = "Error: Function " + arg1 + " returned non-zero error code: " + arg2; break;
	default:	msg = "Error: " + arg1 + " " + arg2 + " " + arg3;
	}
	msgText(msg.c_str());
	return errorCode;
}

//Default text message function sends to stdOut
void msgTextStdOut(const char *msg) {
	std::cout << msg << "\n";
}
//void msgTextStdOut(const std::string &msg) {
//	std::cout << msg << "\n";
//}

//Default update message function sends to stdOut
void msgProgressStdOut(const char *msg, int val) {
	std::cout << "\r" << msg << val;
}

//Old default progress message function sends to stdOut
//void msgProgressStdOut(int val) {
//	std::cout << "\r" << val;
//}

//Default GDAL progress message function sends to stdOut, adapted GDALTermProgress implementation from cpl_progress.cpp
int  msgGDALProgressStdOut(double dfComplete, const char *pszMessage, void *pProgressArg) {
	const int nThisTick = std::min(100, std::max(0, static_cast<int>(dfComplete * 100.0)));
	static int nLastTick = -1;
	if (nThisTick < nLastTick && nLastTick >= 99) nLastTick = -1;
	if (nThisTick <= nLastTick)	return 1;
	while (nThisTick > nLastTick) {
		++nLastTick;
		msgProgress("GDAL process percent complete: ", nLastTick);
		if (nThisTick == 100) msgText("");
	}
	return 1;
}

//Default string text message function sends to stdOut
void msgStringStdOut(const std::string &str) {
	std::cout << str;
}


//Define message function pointers to point to default functions
msgTextFP msgText = (msgTextFP)msgTextStdOut;
msgProgressFP msgProgress = (msgProgressFP)msgProgressStdOut;
msgGDALProgressFP msgGDALProgress = (msgGDALProgressFP)msgGDALProgressStdOut;

msgStringFP msgString = (msgStringFP)msgStringStdOut;

