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
	vector<uint32_t> headerOffsetTable;
	vector<uint32_t> sequenceOffsetTable;
public:
	Database();
	void headerInformation();
	void printSequenceName (int place);
	vector<uint32_t> getSequenceOffsetTable();
	
};

#endif



