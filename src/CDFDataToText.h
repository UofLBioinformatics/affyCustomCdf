/*
 * CDFDataToTextFile.h
 * Class to ASCII CDF file
 *
 * Author: Ernur Saka
 */

#ifndef SRC_CDFDATATOTEXT_H_
#define SRC_CDFDATATOTEXT_H_

#include "AffxString.h"
#include "Err.h"
#include "Fs.h"

//
#include "Parameter.h"
#include "GenericFileReader.h"
#include "StringUtils.h"
#include "CDFFileData.h"
//
#include <fstream>
//

using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_parameter;




class CDFDataToTextFile
{
public:

	CDFDataToTextFile();
    virtual ~CDFDataToTextFile();

    bool open(const AffxString& strFileName);

    bool isOpen();
    void write(const char *psz);
    void writeLine(const AffxString& str);
    void close(void);

protected:

    std::fstream* m_pfileStream;
    char* m_cb;

};

class CDFDataToText
{
	protected:
		AffxString m_strFileNameOut;
		CDFDataToTextFile m_file;
		bool m_bHeader;
		int numCols;
		int maxUnit;
		affxcdf::CCDFFileData* cdf;
		affxcdf::CCDFFileHeader header;

public:
	CDFDataToText(affxcdf::CCDFFileData* _cdf)
	{
		m_bHeader = false;
		numCols = -1;
		maxUnit = -1;
		cdf = _cdf;
	}

	void CreateCDFFile(const std::string& strFileNameOut,
		   int _maxUnit,
           bool bHeader,
           bool bBody);

protected:
	void error(const AffxString& strMessage);
	/*
		 * Output the QC Units
	*/
	void OutputQCUnits();

	/*
	 * Output the file header parameters.
	 */
	void OutputHeader();

	/*
	 * Output the group information.
	 */
	void OutputGroupInfo(GenericData &gdata, const std::string &groupName, int groupIndex);

	/*
	 * Output the groups and sets.
	 */
	void OutputGroupsAndSets(GenericData &gdata);

};


#endif /* SRC_CDFDATATOTEXT_H_ */
