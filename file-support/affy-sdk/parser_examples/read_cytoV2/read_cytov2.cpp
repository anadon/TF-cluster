
#include "calvin_files/fusion/src/FusionCHPMultiDataData.h"
#include "util/convert.h"
#include <iostream>

using namespace std;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_parameter;
using namespace affymetrix_calvin_data;

static void printParams(std::list<ParameterNameValueType> &params)
{
	for (std::list<ParameterNameValueType>::iterator it=params.begin(); it!=params.end(); it++)
	{
		ParameterNameValueType &p = *it;
		wcout << p.GetName() << "=" << p.ToString() << endl;
	}
}

static void printParams(std::vector<ParameterNameValueType> &params)
{
	for (std::vector<ParameterNameValueType>::iterator it=params.begin(); it!=params.end(); it++)
	{
		ParameterNameValueType &p = *it;
		wcout << p.GetName() << "=" << p.ToString() << endl;
	}
}

int main(int argc, char **argv)
{
	const char* fileName = (argc == 1 ? "c:\\temp\\test.cychp" : argv[1]);
	
	try
	{
		FusionCHPData *chp = FusionCHPDataReg::Read(fileName);
		if (chp == NULL)
		{
			cout << "Error reading the file." << endl;
			return 0;
		}

		// Now test to see if its a multi-data file.
		FusionCHPMultiDataData *mchp = FusionCHPMultiDataData::FromBase(chp);
		if (mchp == NULL)
		{
			cout << "CHP file not compatible with this example." << endl;
			delete chp;
			return 1;
		}

		wcout << "File name = " << fileName << endl;
        wcout << "Array type = " << mchp->GetArrayType() << endl;
        wcout << "Alg version = " << mchp->GetAlgVersion() << endl;
        wcout << "Alg name = " << mchp->GetAlgName() << endl;
		printParams(mchp->GetAlgParams());
		printParams(mchp->GetSummaryParams());

		MultiDataType dataTypes1[] = {
			ChromosomeSummaryMultiDataType
		};
		std::string dataTypes1str[] = {
			"ChromosomeSummaryMultiDataType"
		};

		MultiDataType dataTypes2[] = {
			SegmentCNMultiDataType,
			SegmentLOHMultiDataType,
			SegmentCNNeutralLOHMultiDataType,
			SegmentNormalDiploidMultiDataType,
			SegmentMosaicismMultiDataType,
			SegmentNoCallMultiDataType
		};
		std::string dataTypes2str[] = {
			"SegmentCNMultiDataType",
			"SegmentLOHMultiDataType",
			"SegmentCNNeutralLOHMultiDataType",
			"SegmentNormalDiploidMultiDataType",
			"SegmentMosaicismMultiDataType",
			"SegmentNoCallMultiDataType"
		};


        MultiDataType dataTypes3[] = {
			CopyNumberMultiDataType
        };
		std::string dataTypes3str[] = {
			"CopyNumberMultiDataType",
		};

        MultiDataType dataTypes4[] = {
            MarkerABSignalsMultiDataType
        };
        std::string dataTypes4str[] = {
			"MarkerABSignalsMultiDataType",
		};

        MultiDataType dataTypes5[] = {
            AllelePeaksMultiDataType
        };
		std::string dataTypes5str[] = {
			"AllelePeaksMultiDataType",
		};

		 MultiDataType dataTypes6[] = {
            CytoGenotypeCallMultiDataType
        };
		std::string dataTypes6str[] = {
			"CytoGenotypeCallMultiDataType",
		};


		for (int i = 0; i < sizeof(dataTypes1)/sizeof(MultiDataType); i++)
		{
			int n = mchp->GetEntryCount(dataTypes1[i]);
			cout << "#Entries(" << dataTypes1str[i] << ")=" << n << endl;
			ChromosomeMultiDataSummaryData d;
			for (int j = 0; j < min(3,n); j++)
			{
				mchp->GetChromosomeSummaryEntry(dataTypes1[i], j, d);
				cout << "Chr = " << (int)d.chr << endl;
				cout << "Display = " << d.display << endl;
				cout << "StartIndex = " << d.startIndex << endl;
				cout << "MarkerCount = " << d.markerCount << endl;
				cout << "MinSignal = " << d.minSignal << endl;
				cout << "MaxSignal = " << d.maxSignal << endl;
				cout << "MedianCnState = " << d.medianCnState << endl;
				cout << "HetFrequency = " << d.hetFrequency << endl;
				cout << "HomFrequency = " << d.homFrequency << endl;
				printParams(d.metrics);
			}
		}

		for (int i = 0; i < sizeof(dataTypes2)/sizeof(MultiDataType); i++)
		{
			int n = mchp->GetEntryCount(dataTypes2[i]);
			cout << "#Entries(" << dataTypes2str[i] << ")=" << n << endl;
			ChromosomeSegmentData d;
			for (int j = 0; j < min(3,n); j++)
			{
				mchp->GetChromosomeSegmentEntry(dataTypes2[i], j, d);
				cout << "Chr = " << (int)d.chr << endl;
				cout << "SegmentId = " << d.segmentId << endl;
				cout << "MarkerCount = " << d.markerCount << endl;
				cout << "StartPosition = " << d.startPosition << endl;
				cout << "StopPosition = " << d.stopPosition << endl;
				cout << "MeanMarkerDistance = " << d.meanMarkerDistance << endl;
				printParams(d.metrics);
			}
		}

        for (int i = 0; i < sizeof(dataTypes3)/sizeof(MultiDataType); i++) {
            int n = mchp->GetEntryCount(dataTypes3[i]);
            cout << "#Entries(" << dataTypes3str[i] << ")=" << n << endl;
            ProbeSetMultiDataCopyNumberData d;
            for (int j = 0; j < min(n, 3); j++) {
                mchp->GetCopyNumberEntry(dataTypes3[i], j, d);
                cout << "Name = " << d.name << endl;
				cout << "Chr = " << (int)d.chr << endl;
				cout << "Position = " << d.position << endl;
                printParams(d.metrics);
            }
        }
		
		MarkerABSignals d;
        for (int i = 0; i < sizeof(dataTypes4)/sizeof(MultiDataType); i++) {
            int n = mchp->GetEntryCount(dataTypes4[i]);
            cout << "#Entries(" << dataTypes4str[i] << ")=" << n << endl;
           
            for (int j = 0; j < min(n, 3); j++) {
                mchp->GetMarkerABSignalsEntry(dataTypes4[i], j, d);
                cout << "Index = " << d.index << endl;
				printParams(d.metrics);
            }
        }

        for (int i = 0; i < sizeof(dataTypes5)/sizeof(MultiDataType); i++) {
            int n = mchp->GetEntryCount(dataTypes5[i]);
            cout << "#Entries(" << dataTypes5str[i] << ")=" << n << endl;
            AllelePeaks d;
            for (int j = 0; j < min(n, 3); j++) {
				mchp->GetAllelePeakEntry(dataTypes5[i], j, d);
                cout << "Name = " << d.name << endl;
                cout << "Chr = " << (int)d.chr << endl;
				cout << "Position = " << d.position << endl;
				printParams(d.peaks);
            }
        }

		for (int i = 0; i < sizeof(dataTypes6)/sizeof(MultiDataType); i++) {
            int n = mchp->GetEntryCount(dataTypes6[i]);
            cout << "#Entries(" << dataTypes6str[i] << ")=" << n << endl;
            CytoGenotypeCallData d;
            for (int j = 0; j < min(n, 9); j++) {
				mchp->GetCytoGenotypeEntry(dataTypes6[i], j, d);
                cout << "Index = " << d.index << endl;
                cout << "Call = " << (int)d.call << endl;
				cout << "Confidence = " << d.confidence << endl;
				string f = ToStr(d.confidence);
				cout <<  "strConfidence = " << f << endl;
				cout << "ej1 = " << .000001f << endl;
				cout << "ej2 = " << .0000001f << endl;
				cout << "ForcedCall = " << (int)d.forcedCall << endl;
				cout << "ASignal = " << d.aSignal << endl;
				cout << "BSignal = " << d.bSignal << endl;
				cout << "SignalStrength = " << d.signalStrength << endl;
				cout << "Contrast = " << d.contrast << endl;
				printParams(d.metrics);
            }
        }



		delete chp;
	}
	catch (...)
	{
		cout << "Error reading the file." << endl;
	}
	return 0;
}
