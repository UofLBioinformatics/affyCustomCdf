/*
 * CustomCdf.cpp
 *
 * Author: Ernur Saka
 */

#include <iostream>

#include "customCdf.h"
#include "FusionCDFData.h"

#include <map>
#include <cassert>
#include <string>


using namespace std;
using namespace affxcdf;


CustomCdf::CustomCdf(const char * cdfFileName, const char * _reportFileName,
                     unsigned short _controlProbeSetNumber,int _probeLengt,
                     unsigned short _unitType){

	originalCdfFileName = cdfFileName;
  reportFileName = _reportFileName;
	controlProbeSetNumber = _controlProbeSetNumber;
	probeLengt = _probeLengt;
	unitType = _unitType;
	MM = true;
}

CustomCdf::~CustomCdf(){

  PMProbes.clear();
  MMProbes.clear();
  histogramP.clear();

}
/*
 * SetProbeInfo READS PM AND MM PROBES INFO FROM ORIGINAL CDF
 * AND SAVES THEM INTO TWO MAPS IN ORDER PMProbes, MMProbes
 * IT ALSO ASSIGN VALUE TO MM BASED ON THE SIZE OF MMProbes
 */

bool CustomCdf::SetProbeInfo(){

  FusionCDFFileHeader header;
  FusionCDFProbeInformation probeInfo;


  FusionCDFQCProbeSetInformation qcunit;
  FusionCDFQCProbeInformation qcprobe;


  cdf.SetFileName(originalCdfFileName);

  Rprintf("Attempting to read original CDF File: %s\n", cdf.GetFileName().c_str());

  try{

  if(cdf.Read() == false)
  {
    REprintf("Failed to read the CDF file.\n");
    return false;
  }
  }catch (...){
    REprintf("Failed to read cdf file and throw an error by.\n");
    return false;
  }

  int nsets = cdf.GetHeader().GetNumProbeSets();

  nRows = cdf.GetHeader().GetRows();
  nCols = cdf.GetHeader().GetCols();

  for(int iset=0; iset<nsets; iset++) {

    FusionCDFProbeSetInformation set;
    cdf.GetProbeSetInformation(iset, set);
    int ngroups = set.GetNumGroups();

    for (int igroup=0; igroup<ngroups; igroup++)
    {
      FusionCDFProbeGroupInformation group;
      set.GetGroupInformation(igroup, group);
      int ncells = group.GetNumCells();

      bool _MM = false;
      if(group.GetNumCellsPerList() == 2)
        _MM = true;

      for (int icell=0; icell<ncells; icell++)
      {
        probe pInfo;
        group.GetCell(icell,probeInfo);

        pInfo.tbase = probeInfo.GetTBase();
        pInfo.pbase = probeInfo.GetPBase();
        pInfo.x = probeInfo.GetX();
        pInfo.y = probeInfo.GetY();

        int hashKey = (pInfo.x + 1) + (nRows*pInfo.y);
        PMProbes.insert(std::pair<int, probe>(hashKey,pInfo));

        if(_MM == true){

          group.GetCell(icell+1,probeInfo);
          pInfo.tbase = probeInfo.GetTBase();
          pInfo.pbase = probeInfo.GetPBase();
          pInfo.x = probeInfo.GetX();
          pInfo.y = probeInfo.GetY();

          MMProbes.insert(std::pair<int, probe>(hashKey,pInfo));
          icell++;
        }
      }
    }
  }

  if(MMProbes.size() != 0)
    MM = true;
  else
    MM = false;

  return true;
}

/*
 * CHECK SIZE OF PMProbes AND MMProbes
 */

bool CustomCdf::CheckSizeOfProbeMaps (){

  if(PMProbes.size() <= 0) {
    REprintf("PM probes info could not be read from original CDF\n");
    return false;
  }

  if(MM && (PMProbes.size() != MMProbes.size())) {
    REprintf("PM and MM probe numbers are not equal\n");
    return false;
  }

  return true;

}

/*
 * IT CREATES PROBES FROM INPUST GIVEN BY R CODES AND SAVED THEM INTO
 * probesNEW MAP. IT RETURNS THE SIZE OF newProbeSets
 */
int CustomCdf::SetCustomProbes(std::vector<std::string> OriginalProbeSet,
                                SEXP x, SEXP y){

	if(!((OriginalProbeSet.size() == length(x)) && (length(x) == length(y)))){
	  REprintf("Each vector has to have same length \n");
	  return -1;
	} else {
		unsigned int i,len=length(x);

		for(i=0; i< len; i++){
      //probe probesNew{INTEGER(x)[i],INTEGER(y)[i], 'N', 'N'};
      probe probesNew;
      probesNew.x = INTEGER(x)[i];
      probesNew.y = INTEGER(y)[i];
      probesNew.pbase = 'N';
      probesNew.tbase = 'N';
      newProbeSets[std::string(OriginalProbeSet[i])].push_back(probesNew);
		}

		/*
		//Travel Inside the map
		std::map<std::string, std::vector<probe> >::iterator it;

		for (it = newProbeSets.begin(); it != newProbeSets.end(); it++){
			Rprintf("Probe set name %s\n",it->first.c_str());
			std::vector<probe>::iterator itIn;
			for (itIn = it->second.begin(); itIn != it->second.end(); itIn++){
				Rprintf("\t x=%d y=%d\n",(*itIn).x,(*itIn).y);
			}
		}
		 */

		return (newProbeSets.size());

	}

}

/*
 * IT FILLS THE pbase tbase INFO OF NEW PROBE SETS BASED ON THE ORIGINAL CDF
 */

void CustomCdf::FillCustomProbes(const std::string& missiggProbesFile){

	std::map<std::string, std::vector<probe> >::iterator it;

	ofstream outputFile(missiggProbesFile.c_str());

	if(!outputFile){
	  REprintf("Missing probes file could not be opened!\n");
	}
	else{
	  outputFile << "List of probes located in the probe file but not used in the original CDF\n";
	  outputFile << "Name\t\t\t\t\t\tX\tY\n";
	}

	if(newProbeSets.size() != 0) {
  	for (it = newProbeSets.begin(); it != newProbeSets.end(); it++){

  		std::vector<probe>::iterator itIn;

  		for (itIn = it->second.begin(); itIn != it->second.end();){

  			int searchKey = ((*itIn).x +1) + ((*itIn).y * nRows);
  			std::map<int, probe>::iterator it1;
  			it1 = PMProbes.find(searchKey);

  			if(it1!=PMProbes.end())
  			{
  				probe tempProbe = it1->second;
  				(*itIn).pbase = tempProbe.pbase;
  				(*itIn).tbase = tempProbe.tbase;
  				itIn++;
  			}
  			else
  			 {
  			  if(outputFile){
  			    outputFile << it->first.c_str() << "\t" << (*itIn).x << "\t" << (*itIn).y << "\n";
  			  }
  			  if(it->second.size() == 1){
  			     		newProbeSets.erase(it);
  			     		break;
  			  } else
  			    itIn = it->second.erase(itIn);
  			}
  			}

  		}

  	}


	 if(!outputFile)
	   outputFile.close();

	/*
	//Travel Inside the map
	for (it = newProbeSets.begin(); it != newProbeSets.end(); it++){
		Rprintf("Probe set name %s\n",it->first);
		std::vector<probe>::iterator itIn;
		for (itIn = it->second.begin(); itIn != it->second.end(); itIn++){
			Rprintf("\t x=%d y=%d \t pbase=%c tbase=%c\n",(*itIn).x,(*itIn).y,(*itIn).pbase,(*itIn).tbase);
		}
	}
	 */

}

/*
 * REMOVES PROBE SETS WHICH HAS LESS PROBES THAN GIVEN MINIMUM LIMIT
 */
void CustomCdf::RemoveMinProbeSets(int minProbeSet){

	std::map<std::string, std::vector<probe> >::iterator it;

	for (it = newProbeSets.begin(); it != newProbeSets.end(); it++) {
		std::vector<probe>::iterator itIn;
		if(it->second.size() < minProbeSet){
			newProbeSets.erase(it);
		}
	}

}

int CustomCdf::FillIntoCDFProbeSet(affxcdf::CCDFFileData *newCCDFFileData){

	//////////////////////WRITING HEADER PART///////////////////////////////
	newCCDFFileData->m_Header.m_Version = cdf.GetHeader().GetFormatVersion();

	//If m_Header.m_Version > 0 chip type needs to be added
	//m_Header.m_ChipType = ?

	newCCDFFileData->m_Header.m_Rows = nRows;
	newCCDFFileData->m_Header.m_Cols = nCols;
	newCCDFFileData->m_Header.m_NumProbeSets = newProbeSets.size() + controlProbeSetNumber;
	//newCCDFFileData->m_Header.m_NumProbeSets = controlProbeSetNumber + 1;
	//Max probe set number

	if (newCCDFFileData->m_Header.m_Version > 1){
		newCCDFFileData->m_Header.m_NumQCProbeSets = cdf.GetHeader().GetNumQCProbeSets();
		newCCDFFileData->m_Header.m_Reference = cdf.GetHeader().GetReference();
	}

	////////////////////////WRITING QC PART////////////////////////////////
	// Allocate for the QCProbeSets.
	CCDFQCProbeSetInformation *pQCProbeSet;
	newCCDFFileData->m_QCProbeSets.resize(newCCDFFileData->m_Header.m_NumQCProbeSets);

	int numQCUnitsInFile  = cdf.GetHeader().GetNumQCProbeSets();

	FusionCDFQCProbeSetInformation qcunit;
	FusionCDFQCProbeInformation qcprobe;

	for(int QCi=0; QCi<numQCUnitsInFile; QCi++){

    cdf.GetQCProbeSetInformation((QCi),qcunit);
    int QCcells = qcunit.GetNumCells();

    pQCProbeSet = &newCCDFFileData->m_QCProbeSets[QCi];

    pQCProbeSet->m_QCProbeSetType = qcunit.GetQCProbeSetType();
    pQCProbeSet->m_NumCells = QCcells;

    CCDFQCProbeInformation *pQCCell;
    pQCProbeSet->m_Cells.resize(pQCProbeSet->m_NumCells);

    for(int icell=0; icell<QCcells; icell++){

    	qcunit.GetProbeInformation(icell, qcprobe);

    	pQCCell = &pQCProbeSet->m_Cells[icell];

    	pQCCell->m_X = qcprobe.GetX();;
    	pQCCell->m_Y = qcprobe.GetY();
    	pQCCell->m_PLen = qcprobe.GetPLen();
    	pQCCell->m_Background = 0;
    	pQCCell->m_PMProbe = 0;

    }

	}
  //////////////////////////WRITING PROBE SETS//////////////////////////
  // Allocate for the probe set names.
	int nsets = newProbeSets.size();
  newCCDFFileData->m_ProbeSetNames.Resize(nsets + controlProbeSetNumber);
  newCCDFFileData->m_ProbeSets.resize(nsets + controlProbeSetNumber);

  int iProbeSet=0;
	CCDFProbeSetInformation *pProbeSet = NULL;
	int missMatch = 2;
	bool expectMisMatch = false;

	/////////Write the Control Probe sets from the Original CDF First/////////
	FusionCDFProbeInformation probeInfo;

	for(int iset=0; iset<controlProbeSetNumber; iset++) {

		FusionCDFProbeSetInformation set;
		cdf.GetProbeSetInformation(iset, set);

		pProbeSet = &newCCDFFileData->m_ProbeSets[iProbeSet];

		pProbeSet->m_Index = iProbeSet;
		newCCDFFileData->m_ProbeSetNames.SetName(iProbeSet,
                                           cdf.GetProbeSetName(iset));
		pProbeSet->m_Direction = set.GetDirection();

		pProbeSet->m_NumCellsPerList = set.GetNumCellsPerList();
		pProbeSet->m_NumCells = set.GetNumCells();
		//Probably the unit number check it
		pProbeSet->m_ProbeSetNumber = set.GetProbeSetNumber();
		pProbeSet->m_ProbeSetType = set.GetProbeSetType();
		pProbeSet->m_NumGroups = set.GetNumGroups();

		// Sanity check for relationship of m_NumCellsPerList, m_NumCells and m_NumLists
		if(pProbeSet->m_NumLists != 0 &&
		pProbeSet->m_NumCells / pProbeSet->m_NumLists < 255 &&
		pProbeSet->m_NumCellsPerList != pProbeSet->m_NumCells / pProbeSet->m_NumLists) {
			assert(0 && "CCDFFileData::ReadTextFormat(): m_NumCellsPerList != pProbeSet->m_NumCells / pProbeSet->m_NumLists");
		}

		// If this is an expression probe set and we have 2 cells per list set expectMisMatch flag.
		if(pProbeSet->m_ProbeSetType == ExpressionProbeSetType && pProbeSet->m_NumCellsPerList == 2)
			expectMisMatch = true;

		// Write the blocks.
		CCDFProbeGroupInformation *pBlk;
		pProbeSet->m_Groups.resize(pProbeSet->m_NumGroups);
		for (int iGroup=0; iGroup<pProbeSet->m_NumGroups; iGroup++){

			FusionCDFProbeGroupInformation group;
			set.GetGroupInformation(iGroup, group);
			int ncells = group.GetNumCells();

			pBlk = &pProbeSet->m_Groups[iGroup];
			pBlk->m_GroupIndex = iGroup;
			pBlk->m_ProbeSetIndex = iProbeSet;
			pBlk->m_Name = group.GetName();

			if (pProbeSet->m_ProbeSetType == ExpressionProbeSetType)
			newCCDFFileData->m_ProbeSetNames.SetName(iProbeSet, cdf.GetProbeSetName(iset));

			pBlk->m_NumLists = group.GetNumLists();
			pBlk->m_NumCells = group.GetNumCells();
			pBlk->m_Start = group.GetStart();
			pBlk->m_Stop = group.GetStop();
			pBlk->m_NumCellsPerList = group.GetNumCellsPerList();
			pBlk->m_Direction = group.GetDirection();

			//Start writing cells
			CCDFProbeInformation cell;
			pBlk->m_Cells.resize(ncells);
			unsigned int cellIndex;

			for (int iCell=0; iCell<ncells; iCell++){
				//probe pInfo;
				group.GetCell(iCell,probeInfo);

				cell.m_X = probeInfo.GetX();
				cell.m_Y = probeInfo.GetY();
				cell.m_Expos = cell.m_X + nRows*cell.m_Y;
				cell.m_ProbeLength = probeInfo.GetProbeLength();
				cell.m_PBase = probeInfo.GetPBase();
				cell.m_TBase = probeInfo.GetTBase();
				cell.m_ListIndex = probeInfo.GetListIndex();
				cell.m_ProbeGrouping = probeInfo.GetProbeGrouping();

				if (pProbeSet->m_ProbeSetType == ExpressionProbeSetType)
				{
				cellIndex = (iCell / pProbeSet->m_NumCellsPerList) * pProbeSet->m_NumCellsPerList;
							// If we are expecting pairs of PM/MM then we want the order
							// in m_Cells to be PM first and MM second.
				if (expectMisMatch && cell.m_PBase == cell.m_TBase)
				  ++cellIndex;
				} else {
  				cellIndex = (iCell / pProbeSet->m_NumCellsPerList) * pProbeSet->m_NumCellsPerList;
  				cellIndex += (pProbeSet->m_NumCellsPerList - (iCell % pProbeSet->m_NumCellsPerList) - 1);
				}

				if(cellIndex >= pBlk->m_Cells.size()) {
				assert(0 &&
				 "CCDFFileData::ReadTextFormat(): cellIndex cannot be larger that pBlk->m_Cells.size()" );
				}
				pBlk->m_Cells[cellIndex] = cell;

				if (iCell==0)
				pBlk->m_Start = cell.m_ListIndex;
				else if (iCell == pBlk->m_NumCells-1)
				pBlk->m_Stop = cell.m_ListIndex;

			}//End of for (int iCell=0; iCell<ncells; iCell++)

		}//end of for (int iGroup=0; iGroup<pProbeSet->m_NumGroups; iGroup++)

		iProbeSet++;
	}//end of for(int iset=0; iset<controlProbeSetNumber; iset++)

	int lastProbeSetNumber = 1;
	if(pProbeSet != NULL)
	   lastProbeSetNumber = pProbeSet->GetProbeSetNumber()+1;

	if(MM)
		missMatch = 2;
	else
		missMatch = 1;
	expectMisMatch = false;

	std::map<std::string, std::vector<probe> >::iterator it;
	int numberGroups = 1;

	for (it = newProbeSets.begin(); it != newProbeSets.end(); it++){

		pProbeSet = &newCCDFFileData->m_ProbeSets[iProbeSet];
		pProbeSet->m_Index = iProbeSet;
		newCCDFFileData->m_ProbeSetNames.SetName(iProbeSet, it->first);
		//Get the direction how?????
		std::string pName(it->first);
		if(pName.find("+") != std::string::npos)
			pProbeSet->m_Direction = SenseDirection;
		else if(pName.find("-") != std::string::npos)
			pProbeSet->m_Direction = AntiSenseDirection;
		else
			pProbeSet->m_Direction = NoDirection;


		pProbeSet->m_NumCellsPerList = missMatch;//?????????????????????
		pProbeSet->m_NumCells = it->second.size()*missMatch;

		//Probably the unit number check it
		pProbeSet->m_ProbeSetNumber = lastProbeSetNumber;
		lastProbeSetNumber++;

    typedef enum {
        UNKNOWN_TILE,
        STANDARD_TILE,
        BLOCK_TILE,
        GENE_EXPRESSION_TILE,
        CONTROL_TILE,
        STANDARD_ALTERNATE_TILE,
        STANDARD_VARIANT_TILE,
        UNIVERSAL_TILE,
        COPY_NUMBER_TILE,
        GENOTYPE_CONTROL_TILE,
        EXPRESSION_CONTROL_TILE,
        MARKER_TILE,
        MULTICHANNEL_MARKER_TILE
    } TilingTypes;

	    //int ival = GENE_EXPRESSION_TILE;
	    switch (unitType)
	    {
	    case STANDARD_TILE:
	    case STANDARD_ALTERNATE_TILE:
	    case STANDARD_VARIANT_TILE:
	        pProbeSet->m_ProbeSetType = ResequencingProbeSetType;
	        break;

	    case BLOCK_TILE:
	        pProbeSet->m_ProbeSetType = GenotypingProbeSetType;
	        break;

	    case GENE_EXPRESSION_TILE:
	        pProbeSet->m_ProbeSetType = ExpressionProbeSetType;
	        break;

	    case UNIVERSAL_TILE:
	        pProbeSet->m_ProbeSetType = TagProbeSetType;
	        break;

	    case COPY_NUMBER_TILE:
	        pProbeSet->m_ProbeSetType = CopyNumberProbeSetType;
	        break;

	    case GENOTYPE_CONTROL_TILE:
	        pProbeSet->m_ProbeSetType = GenotypeControlProbeSetType;
	        break;

	    case EXPRESSION_CONTROL_TILE:
	        pProbeSet->m_ProbeSetType = ExpressionControlProbeSetType;
	        break;

	    case MARKER_TILE:
	        pProbeSet->m_ProbeSetType = MarkerProbeSetType;
	        break;

	    case MULTICHANNEL_MARKER_TILE:
	        pProbeSet->m_ProbeSetType = MultichannelMarkerProbeSetType;
	        break;

	    default:
	        pProbeSet->m_ProbeSetType = UnknownProbeSetType;
	        break;
	    }
	    pProbeSet->m_NumGroups = numberGroups;

	    // Sanity check for relationship of m_NumCellsPerList, m_NumCells and m_NumLists
	    if(pProbeSet->m_NumLists != 0 &&
	       pProbeSet->m_NumCells / pProbeSet->m_NumLists < 255 &&
	       pProbeSet->m_NumCellsPerList != pProbeSet->m_NumCells / pProbeSet->m_NumLists) {
	        assert(0 && "CCDFFileData::ReadTextFormat(): m_NumCellsPerList != pProbeSet->m_NumCells / pProbeSet->m_NumLists");
	    }
	    // If this is an expression probe set and we have 2 cells per list set expectMisMatch flag.
	    if(pProbeSet->m_ProbeSetType == ExpressionProbeSetType &&
	       pProbeSet->m_NumCellsPerList == 2)
	        expectMisMatch = true;

		// Write the blocks.
		CCDFProbeGroupInformation *pBlk;
		pProbeSet->m_Groups.resize(pProbeSet->m_NumGroups);

		for (int iGroup=0; iGroup<pProbeSet->m_NumGroups; iGroup++)
		{
			pBlk = &pProbeSet->m_Groups[iGroup];
			pBlk->m_GroupIndex = iGroup;
			pBlk->m_ProbeSetIndex = iProbeSet;

			pBlk->m_Name = it->first;

			 if (pProbeSet->m_ProbeSetType == ExpressionProbeSetType)
				 newCCDFFileData->m_ProbeSetNames.SetName(iProbeSet, it->first);

			pBlk->m_NumLists = it->second.size();
			pBlk->m_NumCells = pBlk->m_NumLists * missMatch;
			pBlk->m_Start = 0;
			pBlk->m_Stop = pBlk->m_NumLists- 1;
			pBlk->m_NumCellsPerList = pProbeSet->m_NumCellsPerList;
			pBlk->m_Direction = pProbeSet->m_Direction;

			//Start writing cells
			CCDFProbeInformation cell;
			pBlk->m_Cells.resize(pBlk->m_NumCells);
			unsigned int cellIndex=0;

			std::vector<probe>::iterator itIn;
			int iCell = 0;
			for (itIn = it->second.begin(); itIn != it->second.end(); itIn++){
			  //PM Probe writing
			  cell.m_X = (*itIn).x;
				cell.m_Y = (*itIn).y;
				cell.m_Expos = cell.m_X + nRows*cell.m_Y;
			  cell.m_ProbeLength = probeLengt;
				cell.m_PBase = (*itIn).pbase;
				cell.m_TBase = (*itIn).tbase;
				cell.m_ListIndex = cellIndex;
        if (pProbeSet->m_ProbeSetType == ExpressionProbeSetType)
        {
          //cellIndex = (iCell / pProbeSet->m_NumCellsPerList) * pProbeSet->m_NumCellsPerList;
                          // If we are expecting pairs of PM/MM then we want the order
                          // in m_Cells to be PM first and MM second.
          if (expectMisMatch && cell.m_PBase == cell.m_TBase)
            ++cellIndex;
          else if(!expectMisMatch)
            ++cellIndex;
        }
        else
        {
          cellIndex = (iCell / pProbeSet->m_NumCellsPerList) * pProbeSet->m_NumCellsPerList;
          cellIndex += (pProbeSet->m_NumCellsPerList - (iCell % pProbeSet->m_NumCellsPerList) - 1);
        }

        if(iCell >= pBlk->m_Cells.size()) {
          assert(0 &&
                 "CCDFFileData::ReadTextFormat(): cellIndex cannot be larger that pBlk->m_Cells.size()" );
        }
        pBlk->m_Cells[iCell] = cell;
        if (iCell==0)
            pBlk->m_Start = cell.m_ListIndex;
        else if (iCell == pBlk->m_NumCells-1)
            pBlk->m_Stop = cell.m_ListIndex;

				iCell++;
				//MM probe writing
				if(MM)
				{

					int searchKey = ((*itIn).x +1) + ((*itIn).y * nRows);
					std::map<int, probe>::iterator it1;
					it1 = MMProbes.find(searchKey);
					if(it1!=PMProbes.end())
					{
						probe tempProbe = it1->second;
						cell.m_X = tempProbe.x;
						cell.m_Y = tempProbe.y;
						cell.m_Expos = cell.m_X + nRows*cell.m_Y;
						cell.m_ProbeLength = probeLengt;
						cell.m_PBase = tempProbe.pbase;
						cell.m_TBase = tempProbe.tbase;
						cell.m_ListIndex = cellIndex;
					}
					else
					  REprintf("Error MM Probe cannot be found %s\n",it->first.c_str());

            if (pProbeSet->m_ProbeSetType == ExpressionProbeSetType)
            {
                //cellIndex = (iCell / pProbeSet->m_NumCellsPerList) * pProbeSet->m_NumCellsPerList;
                                // If we are expecting pairs of PM/MM then we want the order
                                // in m_Cells to be PM first and MM second.
                if (expectMisMatch && cell.m_PBase == cell.m_TBase)
                    ++cellIndex;
            }
            else
            {
                cellIndex = (iCell / pProbeSet->m_NumCellsPerList) * pProbeSet->m_NumCellsPerList;
                cellIndex += (pProbeSet->m_NumCellsPerList - (iCell % pProbeSet->m_NumCellsPerList) - 1);
            }

            if(iCell >= pBlk->m_Cells.size()) {
              assert(0 &&
                     "CCDFFileData::ReadTextFormat(): cellIndex cannot be larger that pBlk->m_Cells.size()" );
            }

            pBlk->m_Cells[iCell] = cell;

            if (iCell==0)
                pBlk->m_Start = cell.m_ListIndex;
            else if (iCell == pBlk->m_NumCells-1)
                pBlk->m_Stop = cell.m_ListIndex;

					iCell++;

				}//end of if(MM)
			}
		}

		iProbeSet++;
	}

	return (lastProbeSetNumber -1);
}

/*
 * CALCULATES THE HISTOGRAM OF NUMBER OF PROBES PER PROBE SET
 * AND SAVED IT INTO A VECTOR NAMED histogramP
 */
void CustomCdf::GetHistogram(){

  std::map<std::string, std::vector<probe> >::iterator it;
  std::map<int, int>::iterator itH;

  for (it = newProbeSets.begin(); it != newProbeSets.end(); it++){
    itH  = histogramP.find(it->second.size());
    if(itH != histogramP.end())
      itH->second++;
    else{
      histogramP[it->second.size()] = 1;
    }
  }

}
/*
 * PRINTS HISTOGRAM INTO REPORT FILE
 */
void CustomCdf::PrintHistogram(){

  std::map<int, int>::iterator itH;
  unsigned short max = 1;
  ofstream outputFile(reportFileName,std::ofstream::app);

  outputFile << "\nTABLE: Distribution of Number of probes per probe set\n";
  outputFile << "\nProbe Number | Number of Probe Set\n";

  for (itH = histogramP.begin(); itH != histogramP.end(); itH++){
    outputFile << "      " ;
    outputFile.width(6); outputFile << std::left  << itH->first << " | ";
    outputFile << "         ";
    outputFile.width(10); outputFile << std::left  << itH->second << "\n";
  }

  outputFile << "\nHISTOGRAM OF THE NUMBER OF PROBES PER PROBE SETS\n\n";
  outputFile << "X: Number of probes\n";
  outputFile << "Y: Number of probe sets\n";
  outputFile << "Each star represent a probe set that contains indicated number of probes at the Y dimension.\n\n";

  for (itH = histogramP.begin(); itH != histogramP.end(); itH++){

    outputFile.width(4); outputFile << std::left << itH->first;
    if(max < itH->second)
      max = itH->second;
    for(unsigned short i=0; i<itH->second; i++){
      outputFile << "*";
    }
    outputFile << "\n";

  }

  outputFile << "    ";
  for(unsigned short i=0; i<max-3; i++)
    outputFile << "-";

  outputFile << "\n";
  outputFile.width(5); outputFile << std::right << "1" ;

  for(unsigned short i=0; i<max-7; i++)
    outputFile << " ";

  outputFile << max << "\n";
  outputFile.close();

}

/*
void CustomCdf::WriteIntoFile(const char * outputCdfName, unsigned short minProbeSet){

	CDFData data(outputCdfName);
	data.SetArrayRows(nRows);
	data.SetArrayCols(nCols);

	cdf.SetFileName(outputCdfName);
	CDFFileWriter* writer = new CDFFileWriter(*(cdf.GetHeader().getCalvinData()));

	u_int8_t unitType,direction, cellsPerAtom;
	u_int32_t atoms, cells;
	CDFProbeSetWriter* probeWriter;

	writer->OpenDataGroup(L"affy probe set 1", 1);
	unitType = 3;
	direction = 1;
	atoms = 2;
	cellsPerAtom = 2;

	cells = atoms*cellsPerAtom;

	// need to figure out from the previous probe set number
	wchar_t probeSetName[100] = L"Probe Set Name";
	u_int32_t i=0;

	probeWriter = writer->CreateProbeSetWriter(probeSetName,unitType,direction,atoms,cells,i,cellsPerAtom);
	probeWriter->WriteHeader();
	probeWriter->Write(10,10,1,3,'C','A');
	probeWriter->Write(10,11,1,3,'C','A');
	probeWriter->Write(11,10,2,4,'C','A');
	probeWriter->Write(11,11,2,4,'C','A');

	probeWriter->Close();
	delete probeWriter;

	writer->CloseDataGroup();

	writer->OpenDataGroup(L"affy probe set 2", 2);
	unitType = 3;
	direction = 1;
	atoms = 2;
	cellsPerAtom = 2;
	cells = atoms*cellsPerAtom;
	probeWriter =
	   writer->CreateProbeSetWriter(L"xda block name 1", unitType,
									direction, atoms, cells, 0,cellsPerAtom);
	probeWriter->WriteHeader();
	probeWriter->Write(12,12,1,3,'C','A');
	probeWriter->Write(12,13,1,3,'C','A');
	probeWriter->Write(13,13,2,4,'C','A');
	probeWriter->Write(13,12,2,4,'C','A');
	probeWriter->Close();
	delete probeWriter;

	unitType = 3;
	direction = 1;
	atoms = 2;
	cellsPerAtom = 2;
	cells = atoms*cellsPerAtom;
	probeWriter =
	writer->CreateProbeSetWriter(L"xda block name 2", unitType,
										direction, atoms, cells, 0,cellsPerAtom);
	probeWriter->WriteHeader();
	probeWriter->Write(12,12,1,3,'C','A');
	probeWriter->Write(12,13,1,3,'C','A');
	probeWriter->Write(13,13,2,4,'C','A');
	probeWriter->Write(13,12,2,4,'C','A');
	probeWriter->Close();
	delete probeWriter;

	writer->CloseDataGroup();

	delete writer;

	//CDFCntrlFileWriter* qcWriter = new CDFCntrlFileWriter(cdf.GetHeader().getCalvinData());

	//Do it later if possible
	//data.SetProbeSetCnt(nProbeSets, Expression);
/*
	//Write QC probes
	CDFCntrlFileWriter* qcWriter = new CDFCntrlFileWriter(data);

	int numQCUnitsInFile  = cdf.GetHeader().GetNumQCProbeSets();

	FusionCDFQCProbeSetInformation qcunit;
	FusionCDFQCProbeInformation qcprobe;

	for(int QCi=1; QCi<=numQCUnitsInFile; QCi++){

		cdf.GetQCProbeSetInformation((QCi-1),qcunit);
		qcWriter->OpenDataGroup(qcunit.GetQCProbeSetType(),numQCUnitsInFile);
	    int QCcells = qcunit.GetNumCells();

	    m_file.writeLine("CellHeader=X\tY\tPROBE\tPLEN\tATOM\tINDEX");
	    for(int cell=0; cell<QCcells; cell++){

	    	GetCntrlProbeSetWriter* qcProbeWriter;
	    	qcWriter->

	    		qcunit.GetProbeInformation(cell, qcprobe);
	    		int index = (qcprobe.GetY() * numCols) + qcprobe.GetX();
	    		m_file.writeLine("Cell" + ::getInt((int)(cell+1)) + "\t" + ::getInt((int)qcprobe.GetX()) + "\t" \
	    				+ ::getInt((int)qcprobe.GetY()) + "\t" + "N" + \
						"\t" + ::getInt((int)qcprobe.GetPLen()) + "\t" + ::getInt((int)0) + "\t" + ::getInt((int)index));
	    	}
	    	m_file.writeLine("");
	    }

*/


	//CDFFileWriter* writer = new CDFFileWriter(cdf.FusionCDFData);
	/*

	u_int8_t direction, cellsPerAtom;
	u_int32_t atoms, cells;
	CDFProbeSetWriter* probeWriter;

	u_int8_t unitType, direction, cellsPerAtom;
	u_int32_t atoms, cells;
	CDFProbeSetWriter* probeWriter;
		//  const wstring probeSetName, groupSetName;

		    writer->OpenDataGroup(L"affy probe set 1", 1);
		    unitType = 3;
		    direction = 1;
		    atoms = 2;
		    cellsPerAtom = 2;
		    //cells = atoms*cellsPerAtom;

		    cells = atoms*cellsPerAtom;

		    // need to figure out from the previous probe set number
		    wchar_t probeSetName[100] = L"Probe Set Name";
		    u_int32_t i=0;

		    probeWriter = writer->CreateProbeSetWriter(probeSetName,unitType,direction,atoms,cells,i,cellsPerAtom);
		    probeWriter->WriteHeader();
		    probeWriter->Write(10,10,1,3,'C','A');
		    probeWriter->Write(10,11,1,3,'C','A');
		    probeWriter->Write(11,10,2,4,'C','A');
		    probeWriter->Write(11,11,2,4,'C','A');

		    probeWriter->Close();
		    delete probeWriter;

	 */
//}
void CustomCdf::WriteTextCDFFile(affxcdf::CCDFFileData *newCCDFFileData, const std::string& strFileNameOut, int _maxUnit)
{

	CDFDataToText newCdfFile(newCCDFFileData);
	newCdfFile.CreateCDFFile(strFileNameOut,_maxUnit,true,true);

}

void CustomCdf::SetUnitType(unsigned short _unitType){
	unitType = _unitType;
}
void CustomCdf::SetControlProbeSetNumber(unsigned short _controlProbeSetNumber){
	controlProbeSetNumber = _controlProbeSetNumber;
}
/*
void CustomCdf::SetMaxControlProbeSet(unsigned short _maxControlProbeSet){
	maxControlProbeSet = _maxControlProbeSet;
}*/


