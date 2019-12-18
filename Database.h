#ifndef DATABASE_H
#define DATABASE_H

#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

class Database
{
private:

	string name;
	
	vector<uint32_t> headerOffsetTable;
	vector<uint32_t> sequenceOffsetTable;
	int numberOfSequences;
	
	ofstream * pfichier; 
	
public:

	Database(string name, ofstream * pfichier); 
	string getName() const;
	vector<uint32_t> getSequenceOffsetTable() const;
	vector<uint32_t> getHeaderOffsetTable() const;
	int getNumberOfSequences() const;
	
	void indexInformation();
	void printSequenceName (int * places, int * scores);
};

#endif
