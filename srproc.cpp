/**
* SRProcess.cpp : A utility program that PROCesses next-gen multiplexed 
* Sequencing Reads that have been aligned to reference sequences
* and reports read counts per reference sequences
*
* Sequencing read data has the following fields:
*  1-8:		SequencerName	flowCell	lane	tile	x	y	tag	read
*  9-15:	sequence	quality	ref	pos	direction	map	mapQuality
*/

#include "SequencingReadProcessor.h"

void ShowUsage(char* progName);

int main(int argc, char* argv[])
{
	if (argc == 1
		|| argc == 2 && argv[1][0] == '-' && tolower(argv[1][1]) == 'h') {
		ShowUsage(argv[0]);
		return 0;
	}
	else if (argc != 3) {
		cout << "Wrong number of arguments for " << argv[0] << endl << endl;
		ShowUsage(argv[0]);
		return 1;
	}

	bool success = false;

	try {
		SequencingReadProcessor seqReadProcessor;

		success = seqReadProcessor.ProcessSequencingRead(argv[1], argv[2]);

		seqReadProcessor.ReportSequencingReads();
	}
	catch (const std::exception& e) {
		cout << "Error occurred running '" << argv[0] << "'" << endl;
		cout << e.what() << endl;
	}
	catch (...) {
		cout << "Unknown error occurred running '" << argv[0] << "'" << endl;
	}

	if (!success)
		return 1;
	else
		cout << "\nDone processing sequencing reads in " << argv[1] << endl;

	return 0;
}

/** Shows how to use the program
*/
void ShowUsage(char* progName)
{
	string name = progName;
	int index = (int)name.rfind('\\');
	if (index >= 0)
		name = name.substr(index + 1);
	cout << name << endl;
	cout << "PROCessing next generation short Sequencing Reads and " << endl;
	cout << "    count aligned reads for each reference sequences";
	cout << endl << endl;
	cout << "Usage: " << name << " <dirName> <sampleBarcodepFile>";
	cout << endl << endl;
	cout << "    <dirName>  The directory containing sequencing read files";
	cout << endl;
	cout << "    <sampleBarcodepFile>  The sample to barcode mapping file";
	cout << endl << endl;
	cout << "    Output files are created in the same directory";
	cout << endl << endl;

	cout << "Examples:" << endl;
	cout << endl;
	cout << "    " << name << endl;
	cout << "        Print this usage" << endl << endl;
	cout << "    " << name << " run1 run1_sampleBarcode.txt" << endl;
	cout << "        Process sequence reads in directory 'run1'" << endl;
	cout << "        using sample mapping from 'run1_sampleBarcode.txt'";
	cout << "        Reads are reported in 'run1_seqReadCounts.csv'";
	cout << endl << endl;
}
