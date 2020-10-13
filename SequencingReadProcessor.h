#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <dirent.h>
#include <vector>
#include <set>
#include <map>

using namespace std;

typedef vector<int>			IntVec;
typedef vector<IntVec>		IntVec2D;
typedef vector<IntVec2D>		IntVec3D;

// Number of fields in sequencing data
const int NUM_DATA_FIELDS = 15;

// Number of reference sequences reads aligned to
// For Human genome 22 + X + Y + M
const int NUM_REF_SEQUENCES = 25;

// Number of reference sequences reads aligned to
const int NUM_LANES = 8;

// Maximum errors allowed in multiplex tag sequence match
const int MAX_BARCODE_ERROR = 2;

// Number of multiplex sample barcode (
const int NUM_BARCODE_SEQUENCES = 12;

// Illumina's 12 multiplex barcode sequences (1-12)
const string BARCODE_SEQUENCES[] = { "ATCACG", "CGATGT", "TTAGGC", "TGACCA",
								     "ACAGTG", "GCCAAT", "CAGATC", "ACTTGA",
								     "GATCAG", "TAGCTT", "GGCTAC", "CTTGTA" };

// Used to store sequencing read data
struct SequencingData
{
	int lane;
	string barcodeSequence;
	string sequence;
	string quality;
	string refSeq;
	string position;
	string direction;
	string map;
	int mapQ;
};

// Used to store sample to barcode mapping information
struct SampleInfo
{
	int lane;
	int plexIndex;
	string name;
	string barcodeSequence;
};

// Used to store processing information for each lane
struct LaneDataInfo
{
	string seqReadFile;
	int lane;
	size_t numReads;
	size_t numBadReads;
	size_t numBadBarcodes;
};

/**
* This is the main class that does all the sequencing read processing
*/
class SequencingReadProcessor
{

public:
	bool ProcessSequencingRead(char* dataDir, char* sampleBarcode);

	void ReportSequencingReads();

private:
	void InitializeProcessor(char* dataDir, char* sampleBarcode);

	void ReadSampleBarcodeMap();

	void SetLane2barcode2sampleInfo(vector<SampleInfo>& sampleInfos);

	void InitializeSequeningReadCounts();

	bool ParseData(char* line, SequencingData& sd);

	vector<string> SplitLine(char* line, const char* sep = "\t");

	string Trim(const string& str, char c);

	void FindInputFiles();

	int GetSamplePlexIndex(const map<string, SampleInfo>& laneSeqIdx2sampleInfo,
		const string& indexSeq, int maxIndexError);

	bool ValidateSampleInfoHeader(char* headerLine);

	int FindBestAlignment(string::const_iterator sIt, string::const_iterator sEnd,
		string::const_iterator tIt, string::const_iterator tEnd, int err,
		int maxIndexError);

	vector<map<string, SampleInfo>> FindLanesWithSampleRead();

private:
	// Sequencing read data directory specified by the 1st command argument
	string sequencingDataDir;

	// Sample to barcode mapping file specified by the 2nd command argument
	string sampleBarcodeMapFile;

	// The list of sequencing read data file in sequencingDataDir
	vector<string> inputDataFiles;

	// Sample to barcode mapping information parsed from sampleBarcodeMapFile
	vector<SampleInfo> sampleInfos;

	// A vector of maps of barcode to SampleInfo in each lane
	vector<map<string, SampleInfo>> lane2barcode2sampleInfo;

	// Map reference sequence name to index
	map<string, int> reference2index;

	// Allocate space for sequencing read counts per reference/lane/sample
	// using vector implementation for fast access
	IntVec3D seqReadCounts;

	// A vector of processing information and statistics for each lane
	vector<LaneDataInfo> laneDataInfo;
};

