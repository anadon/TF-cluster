#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/fusion/src/FusionCDFData.h"
#include "calvin_files/utils/src/StringUtils.h"
#include <iostream>

using namespace std;
using namespace affymetrix_fusion_io;

/*
 * ReadMultipleIntensitiesAndStdDev function demonstrates how
 * to read intensity values from a CEL file into a vector of floats.  
 * 
 * For multi-channel CEL files, this function also demonstrates
 * how to read data from different channels
 */
void ReadMultipleIntensitiesAndStdDev(FusionCELData* cel)
{
	// get the channels in the CEL file
	WStringVector channels = cel->GetChannels();	

	std::vector<float> floats(3);
	float stdDev;
	if (channels.size() > 1)
	{
		//read intensities in a a channel
		cel->SetActiveDataGroup(channels[0]);		
		cel->GetIntensities(0,floats);
		cout << "The values of the first 3 intensities in the first channel are " 
			<< floats[0] << " " << floats[1] << " " << floats[2] << endl;

		cel->SetActiveDataGroup(channels[1]);		
		cel->GetIntensities(0,floats);
		cout << "The values of the first 3 intensities in the second channel are "
			<< floats[0] << " " << floats[1] << " " << floats[2] << endl;

		//switch back to previous channel  
		cel->SetActiveDataGroup(channels[0]);
		stdDev = cel->GetStdv(0);
		cout << "The value of the standard dev is " << stdDev << endl;

		//switch channels 
		cel->SetActiveDataGroup(channels[1]);
		stdDev = cel->GetStdv(0);
		cout << "The value of the standard dev is " << stdDev << endl;	
	}
	else
	{
		cel->GetIntensities(0,floats);
		cout << "The values of the first 3 intensities are " 
			<< floats[0] << " " << floats[1] << " " << floats[2] << endl;

		stdDev = cel->GetStdv(0);
		cout << "The value of the first standard deviation is " << stdDev << endl;	
	}
}

/*
 * ReadSingleChannelCelFile is based on the original sample code for reading
 * CEL files demonstrating backward compatibility when using the Fusion SDK
 * to read single channel CEL files. 
 */
void ReadSingleChannelCelFile(FusionCELData* cel, FusionCDFData* cdf)
{
	cout << "Intensities read by the function ReadSingleChannelCelFile" << endl;
	int n = cel->GetNumCells();
	float sum = 0;
	n = 3;
	for (int i = 0; i < n; i++)
	{
		sum += cel->GetIntensity(i);
	}
	float avg = sum / n;
	cout << "The average intensity is: " << avg << endl;

	int nsets = cdf->GetHeader().GetNumProbeSets();
	std::string name;
	for (int iset=0; iset<nsets; iset++)
	{
		name = cdf->GetProbeSetName(iset);
		sum = 0;
		FusionCDFProbeSetInformation set;
		cdf->GetProbeSetInformation(iset, set);
		int ngroups = set.GetNumGroups();
		for (int igroup=0; igroup<ngroups; igroup++)
		{
			FusionCDFProbeGroupInformation group;
			set.GetGroupInformation(igroup, group);
			int ncells = group.GetNumCells();
			for (int icell=0; icell<ncells; icell++)
			{
				FusionCDFProbeInformation probe;
				group.GetCell(icell, probe);
				sum += cel->GetIntensity(probe.GetX(), probe.GetY());
			}
		}
		avg = sum / set.GetNumCells();
		cout << "The average probe set intensity (" << name << ") is " << avg << endl;
	}
}

/*
 * ReadSingleOrMultiChannelCelFile demonstrates how to read single 
 * and multi-channel CEL and CDF files using the Fusion SDK. 
 *
 * The call to SetActiveDataGroup identifies the channel to obtain data from. 
 * Setting the active datagroup is required when reading multi-channel CEL files.
 * It is optional when reading single channel CEL files.
 */
void ReadSingleOrMultiChannelCelFile(FusionCELData* cel, FusionCDFData* cdf)
{
	cout << "Intensities read by the function ReadSingleOrMultiChannelCelFile" << endl;
	//Get the channels in the CEL file
	WStringVector channels = cel->GetChannels();
	for (int iChannel = 0; iChannel< channels.size(); iChannel++)
	{
		cout << "Reading channel " << iChannel << ": " << StringUtils::ConvertWCSToMBS(channels[iChannel]) << endl;
	
		//Set the active data group (channel) 
		cel->SetActiveDataGroup(channels[iChannel]);
		int n = cel->GetNumCells();
		float sum = 0;		
		for (int i = 0; i < n; i++)
		{
			sum += cel->GetIntensity(i);
		}
		float avg = sum / n;
		cout << "The average intensity is: " << avg << endl;

		int nsets = cdf->GetHeader().GetNumProbeSets();
		std::string name;
	
		for (int iset=0; iset<nsets; iset++)
		{
			name = cdf->GetProbeSetName(iset);
			sum = 0;
			FusionCDFProbeSetInformation set;
			cdf->GetProbeSetInformation(iset, set);
			int ngroups = set.GetNumGroups();
			for (int igroup=0; igroup<ngroups; igroup++)
			{
				FusionCDFProbeGroupInformation group;
				set.GetGroupInformation(igroup, group);
				int ncells = group.GetNumCells();
				for (int icell=0; icell<ncells; icell++)
				{
					FusionCDFProbeInformation probe;
					group.GetCell(icell, probe);
					sum += cel->GetIntensity(probe.GetX(), probe.GetY());
				}
			}
			avg = sum / set.GetNumCells();
			cout << "The average probe set intensity (" << name << ") is " << avg << endl;
		}
	}
}

int main(int argc, char **argv)
{
	const char* celFileName = argv[1];
	const char* cdfFileName = argv[2];
	FusionCELData cel;
	FusionCDFData cdf;

	try
	{
		cel.SetFileName(celFileName);
		if (cel.Read() == false)
		{
			cout << "Failed to read the file." << endl;
			return 0;
		}

		cdf.SetFileName(cdfFileName);
		if (cdf.Read() == false)
		{
			cout << "Failed to read the CDF file." << endl;
			return 0;
		}

		if (cel.IsMultiColor() == true)
		{
			cout << "The CEL file has multiple channels" << endl;
		}
		else
		{
			cout << "The CEL file has one channel" << endl;		
		}

		ReadSingleOrMultiChannelCelFile(&cel, &cdf);
		ReadMultipleIntensitiesAndStdDev(&cel);

		if (cel.IsMultiColor() == false)
			ReadSingleChannelCelFile(&cel, &cdf);


	}
	catch (...)
	{
		cout << "Error in reading the file.";
	}
}


