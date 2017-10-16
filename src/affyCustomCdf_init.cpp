/*
 * affyCustomCdf_init.cpp
 * This file has enterence function to the C++ code for ASCII CDF file creation
 *
 *
 * Author: Ernur Saka
 */
#include "customCdf.h"
#include "CDFFileData.h"
#include "FusionCDFData.h"

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>


#include <wchar.h>
#include <wctype.h>

#include <iostream>
#include <vector>
#include <map>

using namespace std;
using namespace affymetrix_fusion_io;

using namespace affxcdf;


  extern SEXP CreateCustomFile(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                               SEXP, SEXP, SEXP);


  static const R_CallMethodDef CallEntries[] = {
    {"CreateCustomFile", (DL_FUNC) &CreateCustomFile, 11},
    {NULL, NULL, 0}
  };

  extern "C" void R_unload_affyCustomCdf(DllInfo *info) {  // #nocov start
    // Release resources
  }

  extern "C" void R_init_affyCustomCdf(DllInfo *dll)
  {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
  }

  SEXP CreateCustomFile(SEXP _originalCdfFileName, SEXP _newCdfFileName,
                   SEXP _reportFileName, SEXP probeSetNames,SEXP x, SEXP y,
                   SEXP _unitType,SEXP _controlProbeSetNumber,
                   SEXP _probeLength, SEXP _minProbeSet, SEXP _missiggProbesFile){

    PROTECT(probeSetNames = coerceVector(probeSetNames, STRSXP));
    PROTECT(x = coerceVector(x, INTSXP));
    PROTECT(y = coerceVector(y, INTSXP));
    unsigned short minProbeSet = asInteger(_minProbeSet);

    std::vector<std::string> names;
    int sizeNames = LENGTH(probeSetNames);
    names.resize(sizeNames);

    int i;

    for(i=0; i<sizeNames; i++)
      names[i] = CHAR(STRING_ELT(probeSetNames, i));

    const char * originalCdfFileName = (char *) CHAR(STRING_ELT(_originalCdfFileName, 0));
    const char * newCdfFileName = (char *) CHAR(STRING_ELT(_newCdfFileName, 0));
    const char * reportFileName = (char *) CHAR(STRING_ELT(_reportFileName, 0));
    const char * missiggProbesFileT = (char *) CHAR(STRING_ELT(_missiggProbesFile, 0));

    std::string missiggProbesFile(missiggProbesFileT);

    CustomCdf newCdf(originalCdfFileName, reportFileName,
                 asInteger(_controlProbeSetNumber), asInteger(_probeLength),
                 asInteger( _unitType));

    /*
     * reads original cdf file into PM and MM maps inorder to fill
     * pbase, tbase info of custom CDF
     */

    if(!newCdf.SetProbeInfo()) {
      error("Please check your original CDF");
    }


    bool checkPro = newCdf.CheckSizeOfProbeMaps();
    if(!checkPro){
      error("Please check your original CDF");
    }


   /*
    * fills annotations results into newProbeSets map
    */

   int check = newCdf.SetCustomProbes(names,x,y);

    if(check == 0) {
      Rprintf("Please check your probes mapping file and annotation file\n");
      Rprintf("Non of the provided probes were annotated\n");
    }
    else if(check == -1) {
      error("Please check your probes mapping file and annotation file\n");
    }

    /*
     * fills missing info from original CDF which was saved into maps
     */

    newCdf.FillCustomProbes(missiggProbesFile);

    if(minProbeSet > 1)
     newCdf.RemoveMinProbeSets(minProbeSet);

    //figure out histogram of the numebr of probes per probe set
    newCdf.GetHistogram();
    //draw histogram table and histogram into report file
    newCdf.PrintHistogram();

    affxcdf::CCDFFileData *newCCDFFileData = new CCDFFileData();
    int maxUnit = 0;
    if(newCCDFFileData != NULL)
      maxUnit = newCdf.FillIntoCDFProbeSet(newCCDFFileData);
    else {
      error("Out of memory! File data object could not created. \n");
    }

    newCdf.WriteTextCDFFile(newCCDFFileData,newCdfFileName,maxUnit);

    delete newCCDFFileData;


    /*
 //   Rprintf("Entering into CreateFile function\n");
  //  CDFDataToText newCdfFile(newCCDFFileData);
   // newCdfFile.CreateFile(cdfFileName,"DeneQCLastOut.cdf",33041,true,true);
   // Rprintf("Exit from CreateFile function\n");
    */
    UNPROTECT(3);

    return R_NilValue;

  }

/*

int main(){

	//Test the cdf write to see product is binary or ascii
	//writeCustomCdf("TestCdfFile.cdf");

  return 0;
}*/
