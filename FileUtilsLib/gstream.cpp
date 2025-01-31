/*
gstream.cpp - wraps GDAL Dataset in filestream like interface
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



#include "FileUtilsLib.h"
#include "gdal_priv.h"
#include "gdalwarper.h"
#include "gstream.h"

bool gstream::_tifIsRegistered = false;
GDALDriver * gstream::_gdalTifDriver = NULL;
gstream_type gstream::_defaultFileType = GST_Unknown;

std::string gstream::_defaultTifCompression = "LZW";
int gstream::_defaultTifDataType = GDT_Float32;
int gstream::_defaultnBands = 1;
bool gstream::_defaultOverWrite = true;

#ifndef _WIN32
#	define min std::min
#	define max std::max
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//gstream
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

gstream::gstream()
{
	if (!_tifIsRegistered) {
		GDALRegister_GTiff();
		_gdalTifDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
		_tifIsRegistered = true;
	}

	_fileType = GST_Unknown;
	_isOpen = false;
	_hasData = false;
	_hasHeader = false;
	_state = std::ios::goodbit;

	_fileName = "";
	_fltName = "";
	_hdrName = "";
	_prjName = "";

	_gdalDataset = NULL;
	_gdalRasterband = NULL;
	_gdalDataType = GDT_Unknown;
	_sizeofDataType = 0;
	std::memset(_geoTransform, 0, 6 * sizeof(double));
	_blockWidth = 0;
	_blockHeight = 0;

	_nBands = 1;
	_nCols = 0;
	_nRows = 0;
	_nCells = 0;

	_cellWidth = 0.0;
	_cellHeight = 0.0;
	_cellSize = 0.0;
	_xll = 0.0;
	_yll = 0.0;
	_noData = 0.0;

	_prjString = "";

	_resetPos();
}

gstream::gstream(gstream && x)
{
	if (!_tifIsRegistered) {
		GDALRegister_GTiff();
		_gdalTifDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
		_tifIsRegistered = true;
	}

	_overwrite = x._overwrite;

	_fileType = x._fileType;
	_isOpen = x._isOpen;
	_hasData = x._hasData;
	_hasHeader = x._hasHeader;
	_state = x._state;

	_fileName = x._fileName;
	_fltName = x._fltName;
	_hdrName = x._hdrName;
	_prjName = x._prjName;

	_gdalDataset = x._gdalDataset;
	_gdalRasterband = x._gdalRasterband;
	x._gdalDataset = NULL;
	x._gdalRasterband = NULL;

	_gdalDataType = x._gdalDataType;
	_sizeofDataType = x._sizeofDataType;

	std::memset(_geoTransform, 0, 6 * sizeof(double));
	for (size_t i = 0; i < 6; i++) _geoTransform[i] = x._geoTransform[i];

	_blockWidth = x._blockWidth;
	_blockHeight = x._blockHeight;

	_nBands = x._nBands;
	_nCols = x._nCols;
	_nRows = x._nRows;
	_nCells = x._nCells;

	_cellWidth = x._cellWidth;
	_cellHeight = x._cellHeight;
	_cellSize = x._cellSize;
	_xll = x._xll;
	_yll = x._yll;
	_noData = x._noData;

	_band = x._band;
	_filePos = x._filePos;
	_xPos = x._xPos;
	_yPos = x._yPos;
	_toCopy = x._toCopy;
	_copied = x._copied;
	_nXSize = x._nXSize;
	_nYSize = x._nYSize;

	_prjString = x._prjString;
}

gstream::~gstream()
{
	if (_gdalDataset != NULL) GDALClose((GDALDatasetH)_gdalDataset);
}


void gstream::_resetPos() {
	_band = 1;
	_filePos = 0;
	_xPos = 0;
	_yPos = 0;
	_toCopy = 0;
	_copied = 0;
	_nXSize = 0;
	_nYSize = 0;
}

std::string gstream::fileName() {
	return _fileName;
}

uint gstream::nCols() {
	return _nCols;
}

uint gstream::nRows() {
	return _nRows;
}

ull gstream::nCells() {
	return _nCells;
}

double gstream::cellSize() {
	return _cellSize;
}

double gstream::noData() {
	return _noData;
}

int gstream::sizeofDataType()
{
	return _sizeofDataType;
}

int gstream::dataType()
{
	return int(_gdalDataType);
}

std::string gstream::getProjection() const {
	return _prjString;
}

bool gstream::setProjection(std::string prjString) {
	_prjString = prjString;
	return (_gdalDataset != NULL && _gdalDataset->SetProjection(_prjString.c_str()) == 0);
}

bool gstream::getHeader(uint & nCols, uint & nRows, double & cellSize, float & noData) const {
	nCols = _nCols;
	nRows = _nRows;
	cellSize = _cellSize;
	noData = float(_noData);//TODO GetHeader noData as double
	return true;
}

bool gstream::getHeader(uint & nCols, uint & nRows, double & xll, double & yll, double & cellSize, float & noData) const {
	nCols = _nCols;
	nRows = _nRows;
	xll = _xll;
	yll = _yll;
	cellSize = _cellSize;
	noData = float(_noData);//TODO GetHeader noData as double
	return true;
}

bool gstream::compareHeader(const gstream & other) {
	unsigned int nColsA;
	unsigned int nRowsA;
	double xllA;
	double yllA;
	double cellSizeA;
	float noDataA;
	double xyExtTol = 0.0001;

	if (!other.getHeader(nColsA, nRowsA, xllA, yllA, cellSizeA, noDataA))
		return false;

	return (
		nColsA == _nCols &&
		nRowsA == _nRows &&
		//Added tolerance to comparison for when extent a xor b precision has been rounded
		std::abs(xllA - _xll) < xyExtTol &&
		std::abs(yllA - _yll) < xyExtTol &&
		cellSizeA == _cellSize &&
		noDataA == _noData);
}

bool gstream::computeStatistics()
{
	//TODO Requires testing
	if (_fileType == GST_Flt) {
		//TODO
	}
	else if (_fileType == GST_Tif) {
		int approxOK = 0;
		GDALProgressFunc progress = GDALDummyProgress;
		void *progressData = NULL;
		if (_gdalRasterband == NULL) return false;
		if (_gdalRasterband->ComputeStatistics(approxOK, &(_mins[_band]), &(_maxs[_band]), &(_means[_band]), &(_stdDevs[_band]), progress, progressData) != CPLErr::CE_None) return false;
		if (_gdalRasterband->SetStatistics(_mins[_band], _maxs[_band], _means[_band], _stdDevs[_band]) != CPLErr::CE_None) return false;
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//igstream
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//igstream Constructor
igstream::igstream()
{
	_read = &igstream::_readFail;
	_peek = &igstream::_peekFail;
}

//igstream Constructor
igstream::igstream(std::string filename)
{
	_read = &igstream::_readFail;
	_peek = &igstream::_peekFail;
	open(filename);
}

//igstream Move Constructor
igstream::igstream(igstream && x) : gstream(std::move(x))
{
	_fltStream = std::move(x._fltStream);
	_read = x._read;
	_peek = x._peek;
}

//igstream Destructor
igstream::~igstream()
{
	if (is_open()) close();
}

//Open a flt grid
bool igstream::_openFlt(const std::string & filename)
{
	float noData;
	_fltName = GridFN(filename);
	_hdrName = HeaderFN(filename);
	_prjName = PrjFN(filename);
	_fltStream.open(_fltName, std::ios::in | std::ios::binary);
	_hasData = _fltStream.is_open();
	_hasHeader = ReadGridHeader(_hdrName, _nCols, _nRows, _xll, _yll, _cellSize, noData);
	_state = _fltStream.rdstate();

	_sizeofDataType = sizeof(float);
	std::memset(_geoTransform, 0, 6 * sizeof(double));
	_blockWidth = 0;
	_blockHeight = 0;

	_nBands = 1;
	_nCells = ull(_nCols) * ull(_nRows);
	_cellWidth = _cellSize;
	_cellHeight = _cellSize;
	_noData = double(noData);

	_resetPos();

	return _hasData && _hasHeader && _state == std::ios::goodbit;
}

//Open a tif file
bool igstream::_openTif(const std::string & filename)
{
	_gdalDataset = (GDALDataset *)GDALOpenShared(filename.c_str(), GA_ReadOnly);
	if (_gdalDataset == NULL) return false;
	_gdalRasterband = _gdalDataset->GetRasterBand(1);
	if (_gdalRasterband == NULL) return false;

	_gdalDataType = _gdalRasterband->GetRasterDataType();
	_hasData = _gdalDataset != NULL && _gdalRasterband != NULL && _gdalDataType != GDT_Unknown;
	_sizeofDataType = GDALGetDataTypeSizeBytes(GDALDataType(_gdalDataType));
	_hasHeader = _gdalDataset->GetGeoTransform(_geoTransform) == CPLErr::CE_None;
	_gdalRasterband->GetBlockSize(&_blockWidth, &_blockHeight);

	_nBands = uint(_gdalDataset->GetRasterCount());
	_nCols = uint(_gdalRasterband->GetXSize());
	_nRows = uint(_gdalRasterband->GetYSize());
	_nCells = ull(_nCols) * ull(_nRows);

	_cellWidth = _geoTransform[1];
	_cellHeight = _geoTransform[5];
	_cellSize = _cellWidth; //TODO Not square cells
	_xll = _geoTransform[0];
	_yll = _geoTransform[3] + _nRows * _geoTransform[5];
	_noData = _gdalRasterband->GetNoDataValue();

	_prjString = std::string(_gdalDataset->GetProjectionRef());

	_resetPos();

	return _hasData && _hasHeader;
}


//Open a file as a gstream
void igstream::open(const std::string &filename)
{
	if (_isOpen) close();
	_fileName = filename;
	std::string fileExt = "";
	bool hasExt = GetFileExt(_fileName, fileExt);
	bool fileExists = FileExists(_fileName);
	_fileType = GST_Unknown;

	if (hasExt && fileExt.compare("tif") == 0 && fileExists) _fileType = GST_Tif;
	else if (hasExt && (fileExt.compare("flt") == 0 || fileExt.compare("hdr") == 0) && fileExists) _fileType = GST_Flt;
	else if (FileExists(_fileName + ".tif")) {
		_fileName += ".tif";
		_fileType = GST_Tif;
	}
	else if (FileExists(_fileName + ".flt") && FileExists(_fileName + ".hdr")) _fileType = GST_Flt;
	else {

	}

	if (_fileType != GST_Unknown && _defaultFileType == GST_Unknown)
		_defaultFileType = _fileType;

	//Open a Flt grid
	if (_fileType == GST_Flt) {
		_isOpen = _openFlt(_fileName);
		if (_isOpen) {
			_read = &igstream::_readFlt;
			_peek = &igstream::_peekFlt;
		}
	}
	//Open a tif file
	else if (_fileType == GST_Tif) {
		_isOpen = _openTif(_fileName);
		if (_isOpen) {
			_read = &igstream::_readTif;
			_peek = &igstream::_peekTif;
		}
	}
	else {

	}

	if (_isOpen) _state = goodbit;
	else _state &= failbit;
	return;
}

bool igstream::is_open()
{
	return _isOpen;
}

void igstream::close()
{
	//TODO _state

	_fileName = "";
	if (_fileType == GST_Flt) {
		_fltName = "";
		_hdrName = "";
		_prjName = "";
		if (_fltStream.is_open()) _fltStream.close();
		_state = _fltStream.rdstate();
	}
	else if (_fileType == GST_Tif) {
		GDALClose(GDALDatasetH(_gdalDataset));
		_gdalDataset = NULL;
		_gdalRasterband = NULL;
		_gdalDataType = GDT_Unknown;
	}
	else {

	}
	_state = std::ios::goodbit;

	_read = &igstream::_readFail;
	_peek = &igstream::_peekFail;

	_fileType = GST_Unknown;
	_isOpen = false;
	_hasData = false;
	_hasHeader = false;

	_sizeofDataType = 0;
	std::memset(_geoTransform, 0, 6 * sizeof(double));
	_blockWidth = 0;
	_blockHeight = 0;

	_nBands = 1;
	_nCols = 0;
	_nRows = 0;
	_nCells = 0;
	_cellWidth = 0.0;
	_cellHeight = 0.0;
	_cellSize = 0.0;
	_xll = 0.0;
	_yll = 0.0;
	_noData = 0.0;

	_prjString = "";

	_resetPos();

}

igstream & igstream::_readFlt(char *s, std::streamsize n)
{
	_fltStream.read(s, n);
	_state = _fltStream.rdstate();
	_filePos += _fltStream.gcount() / sizeof(float);
	_xPos = _filePos % ull(_nCols);//'=': conversion from 'ull' to 'uint', possible loss of data
	_yPos = _filePos / ull(_nCols);//'=': conversion from 'ull' to 'uint', possible loss of data
	return (*this);
}

igstream & igstream::_readTif(char *s, std::streamsize n)
{
	_toCopy = n / _sizeofDataType;
	_copied = 0;
	while (_toCopy > 0) {
		if (_xPos == 0 && _toCopy > _nCols) {
			_nXSize = _nCols;
			_nYSize = _toCopy / _nCols;
		}
		else {
			_nXSize = min(_toCopy, _nCols - _xPos);
			_nYSize = 1;
		}
		//TODO check CPLError
		//std::cout <<
		//	"\nxpos\t" << _xPos <<
		//	"\nypos\t" << _yPos <<
		//	"\nxsiz\t" << _nXSize <<
		//	"\nysiz\t" << _nYSize << "\n";

		_gdalRasterband->RasterIO(GF_Read, _xPos, _yPos, _nXSize, _nYSize, &(s[_copied * _sizeofDataType]), _nXSize, _nYSize, GDALDataType(_gdalDataType), _sizeofDataType, _nCols * _sizeofDataType, NULL);
		_copied += _nYSize * _nXSize;
		_toCopy -= _nYSize * _nXSize;
		_yPos += (_nYSize * _nXSize + _xPos) / _nCols;
		_xPos = (_nYSize * _nXSize + _xPos) % _nCols;
		_filePos = _yPos * _nCols + _xPos;
		//_gdalRasterband->AdviseRead(_xPos, _yPos, _nXSize, _nYSize, _nXSize, _nYSize, GDALDataType(_gdalDataType), NULL);//TODO overrun?
	}
	_toCopy = 0;
	_copied = 0;
	return (*this);
}

igstream & igstream::_readFail(char *s, std::streamsize n)
{
	_state &= std::ios::failbit;
	return (*this);
}


igstream & igstream::read(char *s, std::streamsize n)
{
	return (this->*_read)(s, n);
}

igstream & igstream::_peekFlt(char *s, std::streamsize n)
{
	_filePos = _fltStream.tellg() / sizeof(float);
	_fltStream.read(s, n);
	_fltStream.seekg(_filePos * sizeof(float));
	_state = _fltStream.rdstate();
	return (*this);
}

igstream & igstream::_peekTif(char *s, std::streamsize n)
{
	_toCopy = n / _sizeofDataType;//'=': conversion from 'std::streamsize' to 'uint', possible loss of data
	_copied = 0;
	while (_toCopy > 0) {
		if (_xPos == 0 && _toCopy > _nCols) {
			_nXSize = _nCols;
			_nYSize = _toCopy / _nCols;
		}
		else {
			_nXSize = min(_toCopy, _nCols - _xPos);
			_nYSize = 1;
		}
		_gdalRasterband->RasterIO(GF_Read, _xPos, _yPos, _nXSize, _nYSize, &(s[_copied * _sizeofDataType]), _nXSize, _nYSize, GDALDataType(_gdalDataType), _sizeofDataType, _nCols * _sizeofDataType, NULL);
		_copied += _nYSize * _nXSize;
		_toCopy -= _nYSize * _nXSize;
		_yPos += (_nYSize * _nXSize + _xPos) / _nCols;
		_xPos = (_nYSize * _nXSize + _xPos) % _nCols;
		//_gdalRasterband->AdviseRead(_xPos, _yPos, _nXSize, _nYSize, _nXSize, _nYSize, GDALDataType(_gdalDataType), NULL);//TODO overrun?
	}
	_toCopy = 0;
	_copied = 0;
	_yPos = _filePos / _nCols;//'=': conversion from 'ull' to 'uint', possible loss of data
	_xPos = _filePos % _nCols;//'=': conversion from 'ull' to 'uint', possible loss of data
	return (*this);
}

igstream & igstream::_peekFail(char *s, std::streamsize n)
{
	_state &= std::ios::failbit;
	return (*this);
}

igstream & igstream::peek(char *s, std::streamsize n)
{
	return (this->*_peek)(s, n);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//ogstream
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Constructor
ogstream::ogstream()
{
	_headerSet = false;
	_waitingOnHeader = false;
	_gdalCreateOptions = nullptr;
	_gdalCreateOptions = CSLAddNameValue(_gdalCreateOptions, "COMPRESS", _defaultTifCompression.c_str());
	_gdalDataType = _defaultTifDataType;
	_sizeofDataType = GDALGetDataTypeSizeBytes(GDALDataType(_gdalDataType));
	_nBands = _defaultnBands;
	_overwrite = _defaultOverWrite;
	_write = &ogstream::_writeFail;
}

//Constructor initialiser
ogstream::ogstream(std::string filename)
{
	_headerSet = false;
	_waitingOnHeader = false;
	_gdalCreateOptions = nullptr;
	_gdalCreateOptions = CSLAddNameValue(_gdalCreateOptions, "COMPRESS", _defaultTifCompression.c_str());
	_gdalDataType = _defaultTifDataType;
	_sizeofDataType = GDALGetDataTypeSizeBytes(GDALDataType(_gdalDataType));
	_nBands = _defaultnBands;
	_overwrite = _defaultOverWrite;
	_write = &ogstream::_writeFail;
	open(filename);
}

//Move Constructor
//igstream Move Constructor
ogstream::ogstream(ogstream && x) : gstream(std::move(x))
{
	_fltStream = std::move(x._fltStream);
	_headerSet = x._headerSet;
	_waitingOnHeader = x._waitingOnHeader;
	_gdalCreateOptions = x._gdalCreateOptions;
	x._gdalCreateOptions = nullptr;
	_write = x._write;
}

//Destructor
ogstream::~ogstream()
{
	if (is_open()) close();
}

//Open a flt grid
bool ogstream::_openFlt(const std::string & filename)
{
	if (_headerSet) {
		_waitingOnHeader = false;
		_fltName = GridFN(filename);
		_hdrName = HeaderFN(filename);
		_prjName = PrjFN(filename);
		_fltStream.open(_fltName, std::ios::out | std::ios::binary);
		_hasData = _fltStream.is_open();
		_state = _fltStream.rdstate();
		_hasHeader = WriteGridHeader(_hdrName, _nCols, _nRows, _xll, _yll, _cellSize, float(_noData));
		_sizeofDataType = sizeof(float);
		_blockWidth = 0;
		_blockHeight = 0;
		_nBands = 1;
	}
	else _waitingOnHeader = true;

	_resetPos();

	return _hasData && _state == std::ios::goodbit;
}

//Open a tif file
bool ogstream::_openTif(const std::string & filename)
{
	if (_headerSet) {
		_waitingOnHeader = false;
		_gdalDataset = _gdalTifDriver->Create(filename.c_str(), _nCols, _nRows, _nBands, GDALDataType(_gdalDataType), _gdalCreateOptions);
		_gdalRasterband = _gdalDataset->GetRasterBand(1);
		_hasData = _gdalDataset != NULL && _gdalRasterband != NULL && _gdalDataType != GDT_Unknown;
		_hasHeader = _gdalDataset->SetGeoTransform(_geoTransform) == CPLErr::CE_None;
		_gdalRasterband->SetNoDataValue(_noData);
	}
	else _waitingOnHeader = true;

	_resetPos();

	return _hasData && _hasHeader;
}


//Open a file as a gstream
void ogstream::open(const std::string &filename)
{
	if (_isOpen) close();
	_fileName = filename;
	std::string fileExt = "";
	bool hasExt = GetFileExt(_fileName, fileExt);

	if (hasExt) {
		if (fileExt.compare("tif") == 0) _fileType = GST_Tif;
		else if (fileExt.compare("flt") == 0 || fileExt.compare("hdr") == 0) _fileType = GST_Flt;
	}
	else if (_defaultFileType == GST_Tif) {
		_fileType = GST_Tif;
		_fileName += ".tif";
	}
	else if (_defaultFileType == GST_Flt) {
		_fileType = GST_Flt;
		_fileName += ".flt";
	}

	bool canWrite = _overwrite || !FileExists(_fileName);

	//Open a Flt grid
	if (canWrite && _fileType == GST_Flt) {
		_isOpen = _openFlt(_fileName);
		if (_isOpen) {
			_write = &ogstream::_writeFlt;
		}
	}
	//Open a tif file
	//TODO delete existing
	else if (canWrite && _fileType == GST_Tif) {
		_isOpen = _openTif(_fileName);
		if (_isOpen) {
			_write = &ogstream::_writeTif;
		}
	}
	else {

	}

	if (_isOpen) _state = goodbit;
	else _state &= failbit;
	return;
}

//TODO This seems to return false after setHeader or copyHeader. Need to check
bool ogstream::is_open()
{
	return _isOpen;
}

void ogstream::close()
{

	_fileName = "";
	if (_fileType == GST_Flt) {
		_fltName = "";
		_hdrName = "";
		_prjName = "";
		if (_fltStream.is_open()) _fltStream.close();
		_state = _fltStream.rdstate();
	}
	else if (_fileType == GST_Tif) {

		this->computeStatistics();

		GDALClose(GDALDatasetH(_gdalDataset));
		_gdalDataset = NULL;
		_gdalRasterband = NULL;
		_gdalDataType = GDT_Unknown;

		CSLDestroy(_gdalCreateOptions);
	}
	else {

	}
	_write = &ogstream::_writeFail;

	_fileType = GST_Unknown;
	_isOpen = false;
	_hasData = false;
	_hasHeader = false;

	_sizeofDataType = 0;
	std::memset(_geoTransform, 0, 6 * sizeof(double));
	_blockWidth = 0;
	_blockHeight = 0;

	_nBands = 1;
	_nCols = 0;
	_nRows = 0;
	_nCells = 0;
	_cellWidth = 0.0;
	_cellHeight = 0.0;
	_cellSize = 0.0;
	_xll = 0.0;
	_yll = 0.0;
	_noData = 0.0;

	_prjString = "";

	_resetPos();

}

ogstream & ogstream::_writeFlt(const char *s, std::streamsize n)
{
	_fltStream.write(s, n);
	_state = _fltStream.rdstate();
	_filePos += n / sizeof(float);
	_xPos = _filePos % _nCols;//'=': conversion from 'ull' to 'uint', possible loss of data
	_yPos = _filePos / _nCols;//'=': conversion from 'ull' to 'uint', possible loss of data
	return (*this);
}

ogstream & ogstream::_writeTif(const char *s, std::streamsize n)
{
	_toCopy = n / _sizeofDataType;//'=': conversion from 'std::streamsize' to 'uint', possible loss of data
	_copied = 0;
	while (_toCopy > 0) {
		if (_xPos == 0 && _toCopy > _nCols) {
			_nXSize = _nCols;
			_nYSize = _toCopy / _nCols;
		}
		else {
			_nXSize = min(_toCopy, _nCols - _xPos);
			_nYSize = 1;
		}
		_gdalRasterband->RasterIO(GF_Write, _xPos, _yPos, _nXSize, _nYSize, (void *)&(s[_copied * _sizeofDataType]), _nXSize, _nYSize, GDALDataType(_gdalDataType), _sizeofDataType, _nCols * _sizeofDataType, NULL);
		_copied += _nYSize * _nXSize;
		_toCopy -= _nYSize * _nXSize;
		_yPos += (_nYSize * _nXSize + _xPos) / _nCols;
		_xPos = (_nYSize * _nXSize + _xPos) % _nCols;
		_filePos = _yPos * _nCols + _xPos;
	}
	_toCopy = 0;
	_copied = 0;
	return (*this);
}

ogstream & ogstream::_writeFail(const char *s, std::streamsize n)
{
	_state &= std::ios::failbit;
	return (*this);
}


ogstream & ogstream::write(const char *s, std::streamsize n)
{
	return (this->*_write)(s, n);
}

void ogstream::setHeader(const uint & nCols, const uint & nRows, const double & xll, const double & yll, const double & cellSize, const float & noData) {

	_nCols = nCols;
	_nRows = nRows;
	_nCells = ull(_nCols) * ull(_nRows);

	_cellWidth = cellSize;
	_cellHeight = cellSize;
	_cellSize = cellSize;
	_xll = xll;
	_yll = yll;
	_noData = double(noData);

	_geoTransform[0] = _xll;
	_geoTransform[1] = _cellWidth;
	_geoTransform[2] = 0.0;
	_geoTransform[3] = yll + _nRows * _cellHeight;
	_geoTransform[4] = 0.0;
	_geoTransform[5] = -_cellHeight;

	_headerSet = true;
	if (_waitingOnHeader) open(_fileName);
}

//TODO return bool or int for error
void ogstream::copyHeader(const std::string & srcFN) {
	igstream srcGS(srcFN);
	this->copyHeader(srcGS);
}

void ogstream::copyHeader(const gstream & src) {

	float noDataFloat;
	src.getHeader(_nCols, _nRows, _xll, _yll, _cellSize, noDataFloat);

	_nCells = ull(_nCols) * ull(_nRows);
	_cellWidth = _cellSize;
	_cellHeight = _cellSize;
	_noData = double(noDataFloat);

	_geoTransform[0] = _xll;
	_geoTransform[1] = _cellWidth;
	_geoTransform[2] = 0.0;
	_geoTransform[3] = _yll + _nRows * _cellHeight;
	_geoTransform[4] = 0.0;
	_geoTransform[5] = -_cellHeight;

	_headerSet = true;
	if (_waitingOnHeader) open(_fileName);
}

bool ogstream::copyProjection(const gstream & src)
{
	_prjString = src.getProjection();
	return (_gdalDataset != NULL && _gdalDataset->SetProjection(_prjString.c_str()) == 0);
}

void ogstream::setTifOptions(std::string compression, int gdaldataType, int nBands, int blockWidth, int blockHeight)
{
	_gdalCreateOptions = CSLSetNameValue(_gdalCreateOptions, "COMPRESS", compression.c_str());
	_gdalDataType = gdaldataType;
	_sizeofDataType = GDALGetDataTypeSizeBytes(GDALDataType(_gdalDataType));
	_nBands = nBands;
}