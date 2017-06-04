// write_chp.cpp : Defines the entry point for the console application.
//

#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/fusion/src/FusionCELData.h"

// For genotyping files
#include "calvin_files/writers/src/CalvinCHPMultiDataFileWriter.h"

// for expression files
#include "calvin_files/writers/src/CalvinCHPQuantificationFileWriter.h"

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_data;
using namespace affymetrix_fusion_io;

/*
 * This example shows how to write an expression results file with signal estimates (such as RMA or PLIER output)
 * If detection p-values (from DABG) are to be included, use the CHPQuantificationDetection and CHPQuantificationDetectionFileWriter objects.
 */
static void WriteExpressionOutput(const string &celFileName, const string &chpFileName)
{
	// Create the data object specifying the parameters (#results, max probe set name)
	// The following will set the number of probe sets to 20 and the maximum marker name length to 32.
	int numProbeSets = 20; // this is the number of probe sets
	int maxProbeSetName = 32; // don't hard code, must loop over all probe sets to determine this value.

	CHPQuantificationData data(chpFileName);
	data.SetEntryCount(numProbeSets, maxProbeSetName);
	
        // Extract the header of the CEL file and add it to the output CHP file. This is to maintain
	// traceability between input and output file. The header of the CEL file contains the header
	// of the parent scan data file (DAT) and a reference to the sample file (ARR).
	FusionCELData cel;
	wstring chipType;
    	try	
	{
		cel.SetFileName(celFileName.c_str());
        	if(!cel.ReadHeader())
		{
            	// error condition. must handle.
			return;
        	}
        	GenericData *gdata = cel.GetGenericData();
		chipType = cel.GetChipType();
        	if (gdata != NULL)
        	{
            		data.GetFileHeader()->GetGenericDataHdr()->AddParent(*gdata->Header().GetGenericDataHdr());
        	}
        	cel.Close();
    	}
    	catch (...)
    	{
        	// error condition. must handle.
		return;
    	}

	// Set the array type, taken from the CEL file header
	data.SetArrayType(chipType);

	// Add the algorithm name, version and any parameters to the header of the CHP file.
	// Summary metrics (e.g. average signal, control metrics such as 3/5' ratios) typically go in the summary parameters section.
	// The list of summary metrics is dependent on the algorithm.
	ParameterNameValueTypeList paramList;
	ParameterNameValueType param;
	data.SetAlgName(L"The algorithm name");
	data.SetAlgVersion(L"The algorithm version");

	param.SetName(L"AlgThreshold"); // this is fake parameter just as an example
	param.SetValueFloat(1.1f);
	paramList.push_back(param);
	data.AddAlgParams(paramList);
	paramList.clear();

	param.SetName(L"TBD metric"); // this is a fake metric just as an example
	param.SetValueFloat(99.9f);
	paramList.push_back(param);
	data.AddSummaryParams(paramList);
	paramList.clear();

	// Add name value parameters to the header, must include information regarding the application (name, version, company).
	param.SetName(L"program-name");
	param.SetValueText(L"The software name goes here");
    	data.GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
    	param.SetName(L"program-version");
    	param.SetValueText(L"Software version goes here");
    	data.GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
    	param.SetName(L"program-company");
    	param.SetValueText(L"Company name goes here");
    	data.GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
	
	// Create a file writer object
	CHPQuantificationFileWriter *writer = new CHPQuantificationFileWriter(data);
	
	// Go to the start of the data set and write the expression results.
	ProbeSetQuantificationData e;
	writer->SeekToDataSet();
	for (int i=0; i<numProbeSets; i++)
	{
		e.name = "xyz";			// the probe set name
		e.quantification = i+20.0f;	// the signal value
		writer->WriteEntry(e);
	}
	delete writer;
}

static void WriteGenotypeOutput(const string &celFileName, const string &chpFileName)
{	
	// Create the data object specifying the parameters (#marker results, columns, max probe set name)
	// The following will set the number of probe sets (markers) to 20 and the maximum marker name length to 32.
	// It will also set the extra column names to those output by the Birdseed algorithm (used for SNP6 genotyping),
	// these are Signal A, Signal B and Forced Call.
	int numProbeSets = 20;    // this is the number of probe sets (# genotype calls)
	int maxProbeSetName = 32; // don't hard code, must loop over all probe sets to determine this value.
	vector<ColumnInfo> cols;
	FloatColumn aSigColumn(L"Signal A");
	cols.push_back(aSigColumn);
	FloatColumn bSigColumn(L"Signal B");
	cols.push_back(bSigColumn);
	UByteColumn forceColumn(L"Forced Call");
	cols.push_back(forceColumn);
	CHPMultiDataData data(chpFileName);
	data.SetEntryCount(GenotypeMultiDataType, numProbeSets, maxProbeSetName, cols);

    	// Extract the header of the CEL file and add it to the output CHP file. This is to maintain
	// traceability between input and output file. The header of the CEL file contains the header
	// of the parent scan data file (DAT) and a reference to the sample file (ARR).
	FusionCELData cel;
	wstring chipType;
    	try
	{
		cel.SetFileName(celFileName.c_str());
        	if(!cel.ReadHeader())
		{
            		// error condition. must handle.
			return;
        	}
        	GenericData *gdata = cel.GetGenericData();
		chipType = cel.GetChipType();
        	if (gdata != NULL)
        	{
            	data.GetFileHeader()->GetGenericDataHdr()->AddParent(*gdata->Header().GetGenericDataHdr());
        	}
        	cel.Close();
    	}
    	catch (...)
    	{
        	// error condition. must handle.
		return;
    	}

	// Set the array type, taken from the CEL file header
	data.SetArrayType(chipType);

	// Add the algorithm name, version and any parameters to the header of the CHP file.
	// Summary metrics (e.g. %AA calls or call rate) typically go in the summary parameters section.
	// The list of summary metrics is dependent on the algorithm.
	ParameterNameValueTypeList paramList;
	ParameterNameValueType param;
	data.SetAlgName(L"The algorithm name");
	data.SetAlgVersion(L"The algorithm version");

	param.SetName(L"AlgThreshold"); // this is fake parameter just as an example
	param.SetValueFloat(1.1f);
	paramList.push_back(param);
	data.AddAlgParams(paramList);
	paramList.clear();

	param.SetName(L"Call Rate"); // this is a fake metric just as an example
	param.SetValueFloat(99.9f);
	paramList.push_back(param);
	data.AddSummaryParams(paramList);
	paramList.clear();

	// Add name value parameters to the header, must include information regarding the application (name, version, company).
	param.SetName(L"program-name");
    	param.SetValueText(L"The software name goes here");
    	data.GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
    	param.SetName(L"program-version");
    	param.SetValueText(L"Software version goes here");
    	data.GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
    	param.SetName(L"program-company");
    	param.SetValueText(L"Company name goes here");
    	data.GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);

	// Create a file writer object
	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);

	// Initialize the probe set object, allocate space for the extra columns (signal A and B and forced call)
	ParameterNameValueType nv;
	ProbeSetMultiDataGenotypeData gn;
	nv.SetName(L"Signal A");
	nv.SetValueFloat(0.0f);
	gn.metrics.push_back(nv);
	nv.SetName(L"Signal B");
	nv.SetValueFloat(0.0f);
	gn.metrics.push_back(nv);
	nv.SetName(L"Forced Call");
	nv.SetValueUInt8(0);
	gn.metrics.push_back(nv);

	// Go to the start of the data set and write the SNP results.
	writer->SeekToDataSet(GenotypeMultiDataType);
	for (int i=0; i<numProbeSets; i++)
	{
		gn.name = "TBD marker name";
		gn.call = SNP_AA_CALL + i%3;		// the genotype call
		gn.confidence = i + 1.0f;		// the genotype call confidence
		gn.metrics[0].SetValueFloat(i + 2.0f);	// the signal A value
		gn.metrics[1].SetValueFloat(i + 3.0f);	// the signal B value
		gn.metrics[2].SetValueUInt8(SNP_AA_CALL + (i+1)%3); // the forced call
		writer->WriteEntry(gn);
	}
	delete writer;
}

int main(int argc, char* argv[])
{
	string fileType = argv[1];	// the type of CHP file to write
	string celFileName = argv[2];	// the full path to the input CEL file.
	string chpFileName = argv[3];	// the full path to the output CHP file.

	if (fileType == "gt")
		WriteGenotypeOutput(celFileName, chpFileName);
	else if (fileType == "exp")
		WriteExpressionOutput(celFileName, chpFileName);

	return 0;
}

