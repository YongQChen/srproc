#include "SequencingReadProcessor.h"

/**
* This function takes in a directory name (dataDir) and processes all
* sequencing read files in the directory.  
*/
bool SequencingReadProcessor::ProcessSequencingRead(char* dataDir, 
	char* sampleBarcodeFile)
{
	InitializeProcessor(dataDir, sampleBarcodeFile);

	int numFilesProcessed = 0;
	for (auto& inputFile : inputDataFiles)
	{
		string curFile = dataDir;
		curFile += "/";
		curFile += inputFile;

		cout << "Reading " << curFile << " ..." << endl;

		ifstream ifs(curFile.c_str());
		if (!ifs.is_open())
			throw runtime_error(string("Cannot open file ") + curFile);

		char line[1024];

		SequencingData sd;

		size_t numLines = 0;
		size_t numBadLines = 0;
		size_t numUnmatchedIndexSeq = 0;

		while (ifs.getline(line, 1024))
		{
			numLines++;

			if (!ParseData(line, sd)) {
				numBadLines++;
				continue;
			}

			map<string, int>::iterator chrIt = reference2index.find(sd.refSeq);
			int refSeqIndex = (int)reference2index.size();
			if (chrIt == reference2index.end()) 	// Add the reference
				reference2index[sd.refSeq] = refSeqIndex;
			else
				refSeqIndex = chrIt->second;

			int plex = GetSamplePlexIndex(lane2barcode2sampleInfo[sd.lane], 
				sd.barcodeSequence, MAX_BARCODE_ERROR);

			if (plex >= 0)	//good index
				seqReadCounts[refSeqIndex][sd.lane][plex]++;
			else		//Bad index
				numUnmatchedIndexSeq++;
		}

		laneDataInfo[sd.lane] = LaneDataInfo{ inputFile, sd.lane, numLines, 
			numBadLines, numUnmatchedIndexSeq };

		ifs.close();
		numFilesProcessed++;
	}

	cout << "\nNumber of data files processed: " << numFilesProcessed << endl;

	return true;
}

void SequencingReadProcessor::InitializeProcessor(char* dataDir,
	char* sampleBarcode)
{
	sequencingDataDir = dataDir;
	sampleBarcodeMapFile = sampleBarcode;

	FindInputFiles();

	ReadSampleBarcodeMap();

	SetLane2barcode2sampleInfo(sampleInfos);

	// Allocate space for sequencing read counts per reference/lane/sample
	InitializeSequeningReadCounts();

	reference2index.clear();
	laneDataInfo.clear();
	laneDataInfo.resize(NUM_LANES);
}

/**
* Read sample to barcode mapping file with heading 
* Lane #	Sample Name	Lib Conc (pM)	Index #
* 1			"S61,62"		14			"1,2"
* 2			"S62,64,65"		14			"2,4,5"
*/
void SequencingReadProcessor::ReadSampleBarcodeMap()
{
	sampleInfos.clear();

	ifstream ifs(sampleBarcodeMapFile);
	if (!ifs)
		throw runtime_error(
			string("Cannot open sample index configuration file ") 
			+ sampleBarcodeMapFile);

	char line[256];

	// Read header
	if (ifs.getline(line, 255) && !ValidateSampleInfoHeader(line)) 
		throw runtime_error(string("Bad header in ") + sampleBarcodeMapFile);

	// Now read sample names and indices	
	while (ifs.getline(line, 255))
	{
		vector<string> sampleFields = SplitLine(line);
		if (sampleFields.size() < 4)
			continue;
		int laneIndex = atoi(sampleFields[0].c_str()) - 1;
		if (laneIndex < 0 || laneIndex > NUM_LANES) 
			throw runtime_error(string("Bad lane '") + sampleFields[0]
						+ "' specified in " + sampleBarcodeMapFile);

		// Sample names are comma separated, may have quotes around them
		strcpy_s(line, Trim(sampleFields[1], '"').c_str());
		vector<string> sampleNames = SplitLine(line, ",");

		// Barcode indices are comma separated, may have quotes around them
		strcpy_s(line, Trim(sampleFields[3], '"').c_str());
		vector<string> sampleSeqIndices = SplitLine(line, ",");

		for (size_t i = 0; i < sampleNames.size(); i++)
		{
			string sampleName = Trim(sampleNames[i], ' ');
			if (!sampleName.empty() && i < sampleSeqIndices.size()) {
				int barcodeIndex = 
					atoi(Trim(sampleSeqIndices[i], ' ').c_str()) - 1;

				if (barcodeIndex < 0 || barcodeIndex > NUM_BARCODE_SEQUENCES)
					throw runtime_error(string("Bad barcode index '") 
							+ sampleSeqIndices[i]
							+ "' specified in " + sampleBarcodeMapFile);

				SampleInfo si{ laneIndex, 0, sampleName, 
					BARCODE_SEQUENCES[barcodeIndex] };
				sampleInfos.push_back(si);
			}
			else
				break;
		}
	}
}

void SequencingReadProcessor::SetLane2barcode2sampleInfo(
	vector<SampleInfo>& sampleInfos)
{
	lane2barcode2sampleInfo.resize(NUM_LANES);
	for (size_t i = 0; i < sampleInfos.size(); i++) {
		SampleInfo& si = sampleInfos.at(i);

		map<string, SampleInfo>& barcode2sampleInfo 
			= lane2barcode2sampleInfo[si.lane];

		si.plexIndex = (int)barcode2sampleInfo.size();
		barcode2sampleInfo[si.barcodeSequence] = si;
	}
}

/**
* Allocate space for sequence read counts
* It is a 3-dimensional jagged matrix (the last dimention may vary)
*/
void SequencingReadProcessor::InitializeSequeningReadCounts()
{
	seqReadCounts.resize(NUM_REF_SEQUENCES);

	for (IntVec2D& seqCounts : seqReadCounts)
	{
		seqCounts.resize(NUM_LANES);
		for (int i = 0; i < NUM_LANES; i++)
			seqCounts[i].resize(lane2barcode2sampleInfo[i].size());
	}
}

//Validate heading "Lane #\tSample Name\tLib Conc (pM)\tIndex #"
bool SequencingReadProcessor::ValidateSampleInfoHeader(char* headerLine)
{
	vector<string> headers = SplitLine(headerLine);

	if (headers.size() == 4 && headers[0] == "Lane #"
		&& headers[1] == "Sample Name" && headers[3] == "Index #")
		return true;
	else {
		cout << "Bad lane heading in sample index map file. ";
		cout << "Expect 'Lane #\tSample Name\tLib Conc (pM)\tIndex #'" << endl;
		return false;
	}
}

/** Parse sequencing data (line) with the following fields
*	1-8:	SequencerName	flowCell	lane	tile	x	y	tag	read
*	9-15:	sequence	quality	ref	pos	direction	length	mapQuality
* 
* Example:
*	1-8:	HWUSI - EAS000	1	4	103	349	1118	ATCACG	1
*	9:		AAAACCGGAGCTTTTGCTGGGGATATATGCTCCTTC
*	10:		aa`ba_a^`a^a``^`a_]YZX\[X^]W\_ ^ \__]]
*	11-15:	chr10.fa	60486	R	36	37
* 
*	Note: field 7 is 0 if not multiplex
*	Note: field 12 is empty if no contif match
*/
bool SequencingReadProcessor::ParseData(char* line, SequencingData& sd)
{
	vector<string> dataFields = SplitLine(line, "\t");

	if (dataFields.size() != NUM_DATA_FIELDS)
		return false;

	// #3 lane
	sd.lane = atoi(dataFields[2].c_str()) - 1;

	// #7 index
	sd.barcodeSequence = dataFields[6];

	// #9, sequence
	sd.sequence = dataFields[8];

	// #10 quality, 
	sd.quality = dataFields[9];

	// #11 reference sequence
	sd.refSeq = dataFields[10];

	// #12 map position
	sd.position = dataFields[11];

	// #13 direction
	sd.direction = dataFields[12];

	// #14 match
	sd.map = dataFields[13];

	// #15 map score
	sd.mapQ = atoi(dataFields[14].c_str());

	return true;
}

vector<string> SequencingReadProcessor::SplitLine(char* line, const char* sep)
{
	vector<string> dataFields;

	char* next_token = NULL;
	char* tok = strtok_s(line, sep, &next_token);
	int i = 1;

	while (tok != NULL) {
		dataFields.push_back(tok);
		tok = strtok_s(NULL, sep, &next_token);
	}

	return dataFields;
}

// Trim character c from both sides of string str
string SequencingReadProcessor::Trim(const string& str, char c)
{
	size_t first = str.find_first_not_of(c);
	if (string::npos == first)
	{
		return str;
	}
	size_t last = str.find_last_not_of(c);
	return str.substr(first, (last - first + 1));
}

/**
* Finding files in a directory is implemented with Dirent 
* and is compatible for Linux and Microsoft Windows
* [GitHub](https://github.com/tronkko/dirent/releases)
*/
void SequencingReadProcessor::FindInputFiles()
{
	struct dirent* entry;
	DIR* dir = opendir(sequencingDataDir.c_str());

	inputDataFiles.clear();

	if (dir == NULL) 
		throw runtime_error( string("Cannot open dir: ") + sequencingDataDir);

	locale loc;
	while ((entry = readdir(dir)) != NULL) {
		if (entry->d_type == DT_REG) {
			string name = entry->d_name;
			string ext;
			for (auto& s : name.substr(name.length() - 4)) 
				ext += tolower(s);

			if(ext == ".txt")
				inputDataFiles.push_back(name);
		}
	}
	closedir(dir);
}

//Get the sample index from index sequence read
int SequencingReadProcessor::GetSamplePlexIndex(
	const map<string, SampleInfo>& laneSeqIdx2sampleInfo,
	const string& indexSeq, int maxIndexError)
{
	map<string, SampleInfo>::const_iterator it, 
		end = laneSeqIdx2sampleInfo.end();

	if ((it = laneSeqIdx2sampleInfo.find(indexSeq)) != end)
		return it->second.plexIndex;	//perfect match

	//Find best match with index error not exceeding maxIndexError
	int bestIdx = -1;
	int minErr = maxIndexError + 1;
	int numHavingMinErr = 0;

	string::const_iterator sIt = indexSeq.begin(), sEnd = indexSeq.end();
	for (it = laneSeqIdx2sampleInfo.begin(); it != end; it++)
	{
		string::const_iterator tIt = it->first.begin(), tEnd = it->first.end();

		//Recursive search
		int err = FindBestAlignment(sIt, sEnd, tIt, tEnd, 0, maxIndexError);

		if (err <= maxIndexError)
		{
			if (err < minErr)
			{
				minErr = err;
				bestIdx = it->second.plexIndex;
				numHavingMinErr = 1;
			}
			else if (err == minErr)
				numHavingMinErr++;
		}
	}

	if (numHavingMinErr > 1)
		bestIdx = -1;	//two or more best matches, set to unmatched

	return bestIdx;
}
/**
* Find best alignment using recursive search
* Allow slide, insertion and deletion upto the maximum index error
* (maxIndexError)
*/
int SequencingReadProcessor::FindBestAlignment(
	string::const_iterator sIt, string::const_iterator sEnd,
	string::const_iterator tIt, string::const_iterator tEnd, 
	int err, int maxIndexError)
{
	int minErr = maxIndexError + 1;
	while (tIt != tEnd && sIt != sEnd)
	{
		if (*tIt != *sIt)
		{
			err++;
			if (err > maxIndexError)
				break;
			else if ((sIt + 1) != sEnd && *tIt == *(sIt + 1))
			{
				int e = FindBestAlignment(sIt + 1, sEnd, tIt, tEnd, 
					err, maxIndexError);
				if (e == err)		//Rest of the sequences match exactly
					return e;
				else if (e < minErr) //Best so far
					minErr = e;
			}
			else if ((tIt + 1) != tEnd && *(tIt + 1) == *sIt)
			{
				int e = FindBestAlignment(sIt, sEnd, tIt + 1, tEnd, 
					err, maxIndexError);
				if (e == err)		//Rest of the sequences match exactly
					return e;
				else if (e < minErr)	//Best so far
					minErr = e;
			}
		}

		sIt++;
		tIt++;
	}

	if (err < minErr)
		minErr = err;

	return minErr;
}

/**
* Report sequencing input and lane informattion and then
* report sequencing read counts per lane, samples and reference sequences
*/
void SequencingReadProcessor::ReportSequencingReads()
{
	//Find the lanes with sequencing reads to report
	vector<map<string, SampleInfo>> laneWithSampleReads 
		= FindLanesWithSampleRead();

	string outputFile = sequencingDataDir + "_seqReadCounts.csv";
	ofstream ofs(outputFile);

	cout << "\nReport sequencing read to " << outputFile << endl;

	// Report inputs
	ofs << "Processed " << inputDataFiles.size() << " files in directory: ";
	ofs << sequencingDataDir << endl << endl;

	ofs << "Sequencing data information for each Lane" << endl;
	ofs << "=========================================" << endl;

	ofs << "Stats/Lane";
	for (auto& maps : laneWithSampleReads) 
		ofs << "," << "Lane#" << maps.begin()->second.lane + 1;
	ofs << endl;

	ofs << "numReads";
	for (auto& maps : laneWithSampleReads)
		ofs << "," << laneDataInfo[maps.begin()->second.lane].numReads;
	ofs << endl;

	ofs << "numBadReads (No data or incomplete data)";
	for (auto& maps : laneWithSampleReads)
		ofs << "," << laneDataInfo[maps.begin()->second.lane].numBadReads;
	ofs << endl;

	ofs << "numBadBarcode (no match or multi-match)";
	for (auto& maps : laneWithSampleReads)
		ofs << "," << laneDataInfo[maps.begin()->second.lane].numBadBarcodes;
	ofs << endl << endl;

	// Report Counts
	ofs << "Sequencing read counts per Lane/sample for each reference" << endl;
	ofs << "=========================================================" << endl;
		
	// Headers
	stringstream ss;
	ofs << "Lane";
	ss << "Reference/Sample";

	for (auto& maps : laneWithSampleReads) {
		for (auto& p : maps) {
			ofs << "," << "Lane#" << p.second.lane + 1;
			ss << "," << p.first;
		}
	}

	ofs << endl;
	ofs << ss.str() << endl;

	// Read counts
	for (auto& r2i : reference2index) {
		ofs << r2i.first;
		IntVec2D& refSeqReadCounts = seqReadCounts[r2i.second];
		for (auto& maps : laneWithSampleReads) {
			for (auto& s2si : maps) {
				SampleInfo& si = s2si.second;
				ofs << "," << refSeqReadCounts[si.lane][si.plexIndex];
			}
		}
		ofs << endl;
	}

	ofs.close();

	cout << "\nDone reporting " << endl;
}

vector<map<string, SampleInfo>> 
SequencingReadProcessor::FindLanesWithSampleRead()
{
	vector<map<string, SampleInfo>> laneWithSampleReads;
	for (size_t i = 0; i < NUM_LANES; i++) {
		size_t counts = 0;
		for (size_t j = 0; j < reference2index.size(); j++) {
			for (auto plexCounts : seqReadCounts[j][i]) {
				counts += plexCounts;
			}
		}
		if (counts > 0) {
			map<string, SampleInfo> sampleWithReads;
			for (auto const& p : lane2barcode2sampleInfo[i])
				sampleWithReads[p.second.name] = p.second;

			laneWithSampleReads.push_back(sampleWithReads);
		}
	}

	return laneWithSampleReads;
}
