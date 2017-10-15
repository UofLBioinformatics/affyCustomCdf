/*
 * customCdf.h
 * Class to prepare CDF data for custom CDF
 *
 * Author: Ernur Saka
 */

#ifndef CUSTOMCDF_H
#define CUSTOMCDF_H

#include <map>
#include <vector>
#include <string>

#include "FusionCDFData.h"
//#include "CDFFileWriter.h"
#include "CDFFileData.h"
//#include "CDFCntrlFileWriter.h"
#include "CDFFileData.h"
#include "CDFDataToText.h"

#include <R.h>
#include <Rdefines.h>

//#include "FusionCDFData.h"

using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_io;

class CustomCdf{

	protected:
		struct probe
		{
				unsigned short x, y;
				char pbase, tbase;
		};

		FusionCDFData cdf;
		const char* originalCdfFileName;
		const char* reportFileName;
		unsigned short unitType;
		unsigned short controlProbeSetNumber;
		//unsigned short maxControlProbeSet;
		int nRows,nCols;
		int probeLengt;
		bool MM;
		//Keeps the new probe sets
		std::map<std::string, std::vector<probe> > newProbeSets;

		//Keeps the probe info obtained from original CDF
		std::map<int, probe> PMProbes;
		std::map<int, probe> MMProbes;
		std::map<int, int> histogramP;


public:

	CustomCdf(const char * cdfFileName, const char * _reportFileName, unsigned short _controlProbeSetNumber,int _probeLengt, unsigned short _unitType);
  //CustomCdf(const char * cdfFileName, unsigned short _controlProbeSetNumber,int _probeLengt, unsigned short _MM, unsigned short _unitType);
	~CustomCdf();
	//void GetProbeInfo(char* fileName);
	SEXP GetProbeInfo(SEXP fname);
	int SetCustomProbes(std::vector<std::string> OriginalProbeSet,SEXP x, SEXP y);
	bool SetProbeInfo();
	bool CheckSizeOfProbeMaps ();

	void SetUnitType(unsigned short);
	void SetControlProbeSetNumber(unsigned short);
	//void SetMaxControlProbeSet(unsigned short);

	void FillCustomProbes(const std::string&);
	//void WriteIntoFile(const char * outputCdfName, unsigned short minProbeSet);

	int FillIntoCDFProbeSet(affxcdf::CCDFFileData *newCCDFFileData);
	void RemoveMinProbeSets(int minProbeSet);
	void WriteTextCDFFile(affxcdf::CCDFFileData *newCCDFFileData,const std::string& strFileNameOut,int _maxUnit);
	void GetHistogram();
	void PrintHistogram();

};

#endif /* CUSTOMCDF_H*/
