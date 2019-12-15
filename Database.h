#ifndef DATABASE_H
#define DATABASE_H

#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

class Database
{
private:

	string index;
	string header;
	string sequence;
	
	uint32_t numberOfSequences;
	uint32_t version;
	uint32_t dbType;
	uint32_t titleLength;
	string title;
	uint32_t timestampLength;
	string timestamp;
	uint64_t residueCount;
	uint32_t maximumSequence;
	
	vector<uint32_t> headerOffsetTable;
	vector<uint32_t> sequenceOffsetTable;
	
public:
	Database();
	void headerInformation();
	void printSequenceName (int * places, int * scores);
	vector<uint32_t> getSequenceOffsetTable();
	vector<uint32_t> getHeaderOffsetTable();
};

#endif
