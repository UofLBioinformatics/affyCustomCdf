/*
 * BinaryCdfToASCII.cpp
 *
 *
 * Author: Ernur Saka
 */

#include "CDFASCIIFileWriter.h"

#include "ExceptionBase.h"
#include "Fs.h"
#include "Util.h"
#include "FusionCDFData.h"

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_parameter;
using namespace affymetrix_fusion_io;

using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_utilities;



// Constructor.
CDFASCIIFile::CDFASCIIFile()
{
  m_cb = NULL;
  m_pfileStream = NULL;
}

// Destructor.
CDFASCIIFile::~CDFASCIIFile()
{
  if (isOpen()) {
    m_pfileStream->close();
    delete m_pfileStream;
  }
  if (m_cb != NULL) {delete[] m_cb;}
}

bool CDFASCIIFile::isOpen() {
  return ((m_pfileStream != NULL) && (m_pfileStream->is_open()));
}

void CDFASCIIFile::close(void)
{
  if (isOpen()) {
      m_pfileStream->close();
      delete m_pfileStream;
      m_pfileStream =  NULL;
      if (m_cb != NULL) {delete[] m_cb; m_cb = NULL;}
  }
}


// Open the file.
bool CDFASCIIFile::open(const AffxString& strFileName )
{

  m_pfileStream = new std::fstream;

  Fs::aptOpen(*m_pfileStream, strFileName, fstream::out | fstream::binary | fstream::trunc);

  if (!isOpen()) {
    delete m_pfileStream;
    m_pfileStream = NULL;
  }
  return isOpen();
}



void CDFASCIIFile::write(const char *psz) {
  if (m_pfileStream != NULL) {
    m_pfileStream->write(psz, (std::streamsize)strlen(psz));
  }
}

void CDFASCIIFile::writeLine(const AffxString& str)
{
  std::ostringstream line;
  line << str << std::endl;
  write(line.str().c_str());
}



/**
 * Error handling for the class.
 * @param strMessage - The error messagae to report.
 */
void CDFASCIIFileWriter::error(const AffxString& strMessage)
{
	Err::errAbort(strMessage);
}

void CDFASCIIFileWriter::OutputQCUnits(){

	int numQCUnitsInFile  = cdf.GetHeader().GetNumQCProbeSets();
	//numQCUnitsInFile = 2;
	//int numQCUnitsInFile  = 1;

    FusionCDFQCProbeSetInformation qcunit;
    FusionCDFQCProbeInformation qcprobe;

    for(int QCi=1; QCi<=numQCUnitsInFile; QCi++){

    	cdf.GetQCProbeSetInformation((QCi-1),qcunit);
    	int QCcells = qcunit.GetNumCells();

    	m_file.writeLine("[QC" + ::getInt((int)QCi) + "]");
    	m_file.writeLine("Type=" + ::getInt((int)qcunit.GetQCProbeSetType()));
    	m_file.writeLine("NumberCells=" + ::getInt((int)QCcells));
    	m_file.writeLine("CellHeader=X\tY\tPROBE\tPLEN\tATOM\tINDEX");

    	for(int cell=0; cell<QCcells; cell++){

    		qcunit.GetProbeInformation(cell, qcprobe);
    		int index = (qcprobe.GetY() * numCols) + qcprobe.GetX();
    		m_file.writeLine("Cell" + ::getInt((int)(cell+1)) + "=" + ::getInt((int)qcprobe.GetX()) + "\t" \
    				+ ::getInt((int)qcprobe.GetY()) + "\t" + "N" + \
					"\t" + ::getInt((int)qcprobe.GetPLen()) + "\t" + ::getInt((int)0) + "\t" + ::getInt((int)index));
    	}
    	m_file.writeLine("");
    }

}

void CDFASCIIFileWriter::run(FusionCDFData _cdf, const std::string& strFileNameOut, int _maxUnit, bool bHeader, bool bBody)
{
  try 	{
	//AffxString m_strFileNameIn = Fs::convertToUncPath(strFileNameIn);
    m_strFileNameOut = Fs::convertToUncPath(strFileNameOut);

	//cdf.SetFileName(strFileNameIn);

	maxUnit = _maxUnit;

	cdf = _cdf;

    //if (cdf != null) {
      numCols = cdf.GetHeader().GetCols();
      Verbose::out(1, "Writing file: " + m_strFileNameOut);
      if (m_file.open(m_strFileNameOut) ){
        m_bHeader = bHeader;
        affymetrix_calvin_io::GenericData genericData;
        OutputHeader();
        OutputQCUnits();
        //Put something to print QC probes
        if (bBody){
            OutputGroupsAndSets(genericData);
        }
        m_file.close();
     // }
      //else 			{
        //Err::errAbort("A problem occurred while processing file (Cannot open output file) " + m_strFileNameIn);
      //}
    }
    else
    	cout << "Failed to read the CDF file." << endl;
  }

  //When things go wrong see if we can die gracefully here.
  catch(Except &e) {
    Verbose::out(0,"");
    error("Exception caught. "
          "Message is: " + ToStr(e.what()));
  }
  catch(const std::bad_alloc &e) {
    Verbose::out(0,"");
    error("std::bad_alloc caught. "
          "The application has run out of memory, or the calvin file is malformed."
          "Message is: " + ToStr(e.what()));
  }
  catch(affymetrix_calvin_exceptions::CalvinException &ce) {
    Verbose::out(0,"");
    error("affymetrix_calvin_exceptions::CalvinException caught. "
          "Affymetrix GeneChip Command Console library has thrown an exception. "
          "Message is: " + affymetrix_calvin_utilities::StringUtils::ConvertWCSToMBS(ce.Description()));
  }
  catch(const std::exception &e) {
    Verbose::out(0,"");
    error("std::exception caught. "
          "Message is: " + ToStr(e.what()));
  }
  catch(...) {
    Verbose::out(0,"");
    error("Unknown exception caught. "
          "No message is available.");
  }
}

/*
 * Output the file header parameters.
 */
void CDFASCIIFileWriter::OutputHeader()
{
	int gor;

	if (m_bHeader)
	{

			m_file.writeLine("[CDF]" );
			m_file.writeLine("Version=GC" + ::getDouble((double)cdf.GetHeader().GetFormatVersion()) + ".0");
			//gor = cdf.GetHeader().GetFormatVersion();
			m_file.writeLine("\n[Chip]");
			m_file.writeLine("Name=" + m_strFileNameOut);
			m_file.writeLine("Rows=" + ::getInt((int)cdf.GetHeader().GetRows()));
			//gor = cdf.GetHeader().GetRows();
			m_file.writeLine("Cols=" + ::getInt((int)cdf.GetHeader().GetCols()));
			//gor = cdf.GetHeader().GetCols();
			m_file.writeLine("NumberOfUnits=" + ::getInt((int)cdf.GetHeader().GetNumProbeSets()));
			//gor = cdf.GetHeader().GetNumProbeSets();
			m_file.writeLine("MaxUnit=" + ::getInt((int)maxUnit));
			m_file.writeLine("NumQCUnits=" + ::getInt((int)cdf.GetHeader().GetNumQCProbeSets()));
			//gor = cdf.GetHeader().GetNumQCProbeSets();
			m_file.writeLine("ChipReference=\n");
	}
}


/*
 * Output the group information.
 */
void CDFASCIIFileWriter::OutputGroupInfo(GenericData &gdata, const std::string &groupName, int groupIndex)
{

	//Get the probe set information
	FusionCDFProbeSetInformation set;
	cdf.GetProbeSetInformation(groupIndex, set);
	int ngroups = set.GetNumGroups();

	//Parts before blocks of an probe set
	m_file.writeLine("[Unit" + ::getInt((int)(groupIndex+1)) + "]");
	m_file.writeLine("Name=NONE");
	m_file.writeLine("Direction=" + ::getInt((int)set.GetDirection()));
	m_file.writeLine("NumAtoms=" + ::getInt((int)(set.GetNumCells()/set.GetNumCellsPerList())));
	m_file.writeLine("NumCells=" + ::getInt((int)set.GetNumCells()));
	m_file.writeLine("UnitNumber=" + ::getInt((int)(groupIndex+1)));
	int UnitType = -1;
	switch (set.GetProbeSetType()) {

		case affxcdf::ExpressionProbeSetType:
		  UnitType = 3;
		  break;
		case affxcdf::GenotypingProbeSetType:
		  UnitType = 2;
		  break;
		case affxcdf::ResequencingProbeSetType:
		  UnitType = 6;
		  break;
		case affxcdf::TagProbeSetType:
		  UnitType = 7;
		  break;
		case affxcdf::CopyNumberProbeSetType:
		  UnitType = 8;
		  break;
		case affxcdf::GenotypeControlProbeSetType:
		  UnitType = 9;
		  break;
		case affxcdf::ExpressionControlProbeSetType:
		  UnitType = 10;
		  break;
		case affxcdf::MarkerProbeSetType:
		  UnitType = 11;
		  break;
		case affxcdf::MultichannelMarkerProbeSetType:
		  UnitType = 12;
		  break;
		default:
		  UnitType = 0;
		  break;
	}
	m_file.writeLine("UnitType=" + ::getInt((int)UnitType));
	m_file.writeLine("NumberBlocks=" + ::getInt((int)ngroups));
	m_file.writeLine("");


	//For each block
	for (int i=0; i<ngroups; i++)
	{
		FusionCDFProbeGroupInformation group;
		set.GetGroupInformation(i, group);
		int ncells = group.GetNumCells();

		//Write the block info before cells
		m_file.writeLine("[Unit" + ::getInt((int)(groupIndex+1)) + "_Block" + ::getInt((int)(i+1)) + "]");
		m_file.writeLine("Name=" + group.GetName());
		m_file.writeLine("BlockNumber=" + ::getInt((int)(i+1)));
		m_file.writeLine("NumAtoms=" + ::getInt((int)(ncells/group.GetNumCellsPerList())));
		m_file.writeLine("NumCells=" + ::getInt((int)ncells));
		m_file.writeLine("StartPosition=" + ::getInt((int)group.GetStart()));
		m_file.writeLine("StopPosition=" + ::getInt((int)group.GetStop()));

		m_file.write("CellHeader=X\tY\tPROBE\tFEAT\tQUAL\tEXPOS\tPOS\tCBASE\tPBASE\tTBASE\tATOM\tINDEX\tCODONIND\tCODON\tREGIONTYPE\tREGION\n");
		for (int icell=0; icell<ncells; icell++)
		{
			FusionCDFProbeInformation probe;
			group.GetCell(icell,probe);

			int index = probe.GetY() * numCols + probe.GetX();

			m_file.writeLine("Cell" + ::getInt((int)icell) + "=" + ::getInt((int)probe.GetX()) + "\t" + ::getInt((int)probe.GetY()) \
					+ "\tN\t" + "control\t" + group.GetName() + "\t" \
					+ ::getInt((int)probe.GetExpos()) + "\t" + ::getInt((int)13) + "\t" + probe.GetTBase() \
					+ "\t" + probe.GetPBase() + "\t" + probe.GetTBase() + "\t" + \
					::getInt((int)probe.GetListIndex()) + "\t" + ::getInt((int)index) + "\t" + \
					::getInt((int)-1) + "\t" + ::getInt((int)-1) + "\t" + ::getInt((int)99) + "\t");
		}
	}
	m_file.writeLine("");
	m_file.writeLine("");
}

/*
 * Output the groups and sets.
 */
void CDFASCIIFileWriter::OutputGroupsAndSets(GenericData &gdata)
{
	int nsets = cdf.GetHeader().GetNumProbeSets();
	for (int i=0; i< nsets; i++)
	{
		string name = cdf.GetProbeSetName(i);
		OutputGroupInfo(gdata, name, i);
	}
}




