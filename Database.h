#ifndef DATABASE_H
#define DATABASE_H

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <cstring>
#include <sstream>

using namespace std;

class Database
{
private:
	string index;
	string header;
	string sequence;
	uint32_t numberOfSequences;
	uint32_t version;
	vector<uint32_t> headerOffsetTable;
	vector<uint32_t> sequenceOffsetTable;
	//int penalExtension; // R
	//int penalOpen; // Q
public:
	Database();
	void headerInformation();
	void printSequenceName (int place);
	vector<uint32_t> getSequenceOffsetTable();
	vector<uint32_t> getHeaderOffsetTable();
};

#endif
