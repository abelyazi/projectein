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
	
	int * headerOffsetTable;
	int * sequenceOffsetTable;
	int numberOfSequences;
	
	ofstream * pfichier; 
	
public:

	Database(string name, ofstream * pfichier); 
	string getName() const;
	int * getSequenceOffsetTable() const;
	int * getHeaderOffsetTable() const;
	int getNumberOfSequences() const;
	
	void indexInformation();
	void printSequenceName (int * places, int * scores);
};

#endif
