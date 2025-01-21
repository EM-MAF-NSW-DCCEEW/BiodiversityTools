/*
gstream.h - wraps GDAL Dataset in filestream like interface
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

#ifndef GSTREAM_H
#define GSTREAM_H

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <map>
//#include "gdal_priv.h"
//#include "gdalwarper.h"

typedef unsigned int uint;
typedef unsigned long long ull;

class GDALDriver;
class GDALDataset;
class GDALRasterBand;
//typedef enum GDALDataType;

typedef enum {
	GST_Unknown = 0,
	GST_Flt = 1,
	GST_Tif = 2
}gstream_type;

class gstream {
protected:

	static const std::ios_base::iostate goodbit = (std::ios_base::iostate)0x0;
	static const std::ios_base::iostate eofbit = (std::ios_base::iostate)0x1;
	static const std::ios_base::iostate failbit = (std::ios_base::iostate)0x2;
	static const std::ios_base::iostate badbit = (std::ios_base::iostate)0x4;

	//Static members
	static bool _tifIsRegistered;
	static GDALDriver * _gdalTifDriver;
	static gstream_type _defaultFileType;
	static std::string _defaultTifCompression;
	static int _defaultTifDataType;
	static int _defaultnBands;
	static bool _defaultOverWrite;

	bool _overwrite;

	gstream_type _fileType;
	bool _isOpen;
	bool _hasData;
	bool _hasHeader;
	std::ios_base::iostate _state;

	//Flt grid stream
	std::string _fileName;
	std::string _fltName;
	std::string _hdrName;
	std::string _prjName;

	std::string _prjString;

	//Tif grid stream
	GDALDataset * _gdalDataset = NULL;
	GDALRasterBand * _gdalRasterband = NULL;
	int _gdalDataType;
	int _sizeofDataType;
	double _geoTransform[6];
	int _blockWidth;
	int _blockHeight;

	//Dimensions
	uint _nBands;
	uint _nCols;
	uint _nRows;
	ull _nCells;
	double _cellWidth;
	double _cellHeight;
	double _cellSize;
	double _xll;
	double _yll;
	double _noData;

	//Statistics band mapped
	std::map<uint, double> _mins;
	std::map<uint, double> _maxs;
	std::map<uint, double> _means;
	std::map<uint, double> _stdDevs;

	//File IO indexing
	uint _band;
	ull _filePos;
	uint _xPos;
	uint _yPos;
	uint _toCopy;
	uint _copied;
	uint _nXSize;
	uint _nYSize;

	void _resetPos();

	gstream();
	gstream(const gstream &) = delete;
	gstream(gstream && x);
	~gstream();

public:
	std::string fileName();
	uint nCols();
	uint nRows();
	ull nCells();
	double cellSize();
	double noData();
	int sizeofDataType();
	int dataType();
	std::string getProjection() const;
	bool setProjection(std::string prjString);
	bool getHeader(uint & nCols, uint & nRows, double & cellSize, float & noData) const;
	bool getHeader(uint & nCols, uint & nRows, double & xll, double & yll, double & cellSize, float & noData) const;
	bool compareHeader(const gstream & other);

	bool computeStatistics();


};

class igstream : public gstream {
private:
	std::ifstream _fltStream;
	bool _openFlt(const std::string & filename);
	bool _openTif(const std::string & filename);

	igstream &  _readFlt(char *s, std::streamsize n);
	igstream &  _readTif(char *s, std::streamsize n);
	igstream &  _readFail(char *s, std::streamsize n);
	igstream & (igstream::*_read)(char *s, std::streamsize n);

	igstream &  _peekFlt(char *s, std::streamsize n);
	igstream &  _peekTif(char *s, std::streamsize n);
	igstream &  _peekFail(char *s, std::streamsize n);
	igstream & (igstream::*_peek)(char *s, std::streamsize n);

public:
	explicit igstream();
	explicit igstream(std::string filename);
	igstream(const igstream &) = delete;
	igstream(igstream && x);
	~igstream();

	//Open and close
	void open(const std::string & filename);
	bool is_open();
	void close();

	//IO
	igstream & read(char *s, std::streamsize n);
	igstream & peek(char *s, std::streamsize n);

	//TODO
	int readRow(char *s, uint row);
	int readRows(char *s, uint first, uint nRows);
	int readNextRow(char *s);


};

class ogstream : public gstream {
private:
	std::ofstream _fltStream;
	bool _openFlt(const std::string & filename);
	bool _openTif(const std::string & filename);
	ogstream &  _writeFlt(const char *s, std::streamsize n);
	ogstream &  _writeTif(const char *s, std::streamsize n);
	ogstream &  _writeFail(const char *s, std::streamsize n);
	ogstream & (ogstream::*_write)(const char *s, std::streamsize n);

	char ** _gdalCreateOptions;

	bool _headerSet;
	bool _waitingOnHeader;

public:
	explicit ogstream();
	explicit ogstream(std::string filename);
	ogstream(const ogstream & x) = delete;
	ogstream(ogstream && x);
	~ogstream();

	//Open and close
	void open(const std::string & filename);
	bool is_open();
	void close();

	//Header, CS and Tif options
	void setHeader(const uint & nCols, const uint & nRows, const double & xll, const double & yll, const double & cellSize, const float & noData);
	void copyHeader(const std::string & srcFN);
	void copyHeader(const gstream & src);
	void setTifOptions(std::string compression = "LZW", int gdaldataType = 6, int nBands = 1, int blockWidth = 0, int blockHeight = 0);
	bool copyProjection(const gstream & src);

	//IO
	ogstream & write(const char *s, std::streamsize n);

	//TODO
	ogstream & write(const char *s, std::streamsize n, uint band);
	int writeRow(const char *s);
	int writeRow(const char *s, uint row);
	int writeRows(const char *s, uint row, uint nRows);

};

#endif