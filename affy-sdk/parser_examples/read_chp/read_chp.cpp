
#include "calvin_files/fusion/src/FusionCHPData.h"
#include "calvin_files/fusion/src/FusionCHPLegacyData.h"
#include "calvin_files/fusion/src/FusionCHPMultiDataData.h"
#include <iostream>

using namespace std;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_io;

static void ExpressionLegacyExample(FusionCHPLegacyData *legchp)
{
	FusionExpressionProbeSetResults psResults;
	float sum = 0;
	int n = legchp->GetHeader().GetNumProbeSets();
	for (int i = 0; i < n; i++)
	{
		legchp->GetExpressionResults(i, psResults);
		sum += psResults.GetSignal();
	}
	float avg = sum / n;
	cout << "The average signal is: " << avg << endl;
}

static void GenotypeLegacyExample(FusionCHPLegacyData *legchp)
{
}

static void LegacyExample(FusionCHPLegacyData *legchp)
{
    if (legchp->GetHeader().GetAssayType() != FusionExpression)
    {
        ExpressionLegacyExample(legchp);
    }
    else if (legchp->GetHeader().GetAssayType() != FusionGenotyping)
    {
        GenotypeLegacyExample(legchp);
    }
}

static void MultiDataExample(FusionCHPMultiDataData *dchp)
{
    int nsnps = dchp->GetEntryCount(GenotypeMultiDataType);
    cout << "snp count = " << nsnps << endl;
    int aacalls = 0;
    int abcalls = 0;
    int bbcalls = 0;
    int nocalls = 0;

    // Output the # of calls
    for (int i=0; i<nsnps; i++)
    {
        u_int8_t call = dchp->GetGenoCall(GenotypeMultiDataType, i);
        if (call == ALLELE_A_CALL)
            ++aacalls;

        else if (call == ALLELE_B_CALL)
            ++bbcalls;

        else if (call == ALLELE_AB_CALL)
            ++abcalls;

        else
            ++nocalls;
    }
    cout << "AA calls = " << aacalls << endl;
    cout << "AB calls = " << abcalls << endl;
    cout << "BB calls = " << bbcalls << endl;
    cout << "No Calls = " << nocalls << endl;

    // Now output the summary metrics. This should include the call rate and gender call
    ParameterNameValueTypeList params = dchp->GetSummaryParams();
    for (ParameterNameValueTypeList::iterator it=params.begin(); it!=params.end(); it++)
    {
        cout << StringUtils::ConvertWCSToMBS(it->GetName()) << " = " << StringUtils::ConvertWCSToMBS(it->ToString()) << endl;
    }
}

int main(int argc, char **argv)
{
	const char* fileName = argv[1];
	try
	{
		FusionCHPData *chp = FusionCHPDataReg::Read(fileName);
		if (chp == NULL)
		{
			cout << "Error reading the file." << endl;
			return 0;
		}

        // Test to see if this is a legacy file
		FusionCHPLegacyData *legchp = FusionCHPLegacyData::FromBase(chp);
		if (legchp != NULL)
		{
            LegacyExample(legchp);
            delete legchp; // chp is just a pointer to legchp so only delete one.
            return 0;
        }
        
        // Now test to see if its a multi-data file.
        FusionCHPMultiDataData *dchp = FusionCHPMultiDataData::FromBase(chp);
        if (dchp != NULL)
        {
            MultiDataExample(dchp);
            delete dchp; // chp is just a pointer to legchp so only delete one.
            return 0;
        }

        // Not multi data or legacy. The project must have included parsers for other
        // CHP file types.
        cout << "CHP file not compatible with this example." << endl;
        delete chp;
	}
	catch (...)
	{
		cout << "Error reading the file." << endl;
	}
	return 0;
}
