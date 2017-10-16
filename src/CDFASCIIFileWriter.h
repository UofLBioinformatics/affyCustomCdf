/*
 * CDFASCIIFileWriter.h
 *
 *  Created on: Mar 28, 2017
 *      Author: Ernur
 */

#ifndef SRC_CDFASCIIFILEWRITER_H_
#define SRC_CDFASCIIFILEWRITER_H_

#include "AffxString.h"
#include "Guid.h"
#include "Util.h"
//
#include "Parameter.h"
#include "GenericFileReader.h"
#include "AffyStlCollectionTypes.h"
#include "StringUtils.h"

#include "FusionCDFData.h"
//
#include <fstream>
//

using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_parameter;

using namespace affymetrix_fusion_io;



class CDFASCIIFile
{
public:

	CDFASCIIFile();
    virtual ~CDFASCIIFile();

    bool open(const AffxString& strFileName);

    bool isOpen();
    void write(const char *psz);
    void writeLine(const AffxString& str);
    void close(void);

protected:

    std::fstream* m_pfileStream;
    char* m_cb;

};

class CDFASCIIFileWriter
{
	protected:
		AffxString m_strFileNameOut;
		CDFASCIIFile m_file;
		bool m_bHeader;
		int numCols;
		int maxUnit;
		FusionCDFData cdf;
		FusionCDFFileHeader header;

public:
		CDFASCIIFileWriter()
	{
		m_bHeader = false;
		numCols = -1;
		maxUnit = -1;

	}

	void run(FusionCDFData _cdf, const std::string& strFileNameOut,
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



#endif /* SRC_CDFASCIIFILEWRITER_H_ */
