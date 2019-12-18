#ifndef DATABASE_H
#define DATABASE_H

#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

class Database
{
private:

	string indexFileName = "./uniprot_sprot.fasta.pin";
	string sequenceFileName = "./uniprot_sprot.fasta.psq";
	string headerFileName = "./uniprot_sprot.fasta.phr";
	
	int * headerOffsetTable;
	int * sequenceOffsetTable;
	int numberOfSequences;
	
	ofstream * pfichier; 
	
public:

	Database(ofstream * pfichier); 
	
	string getIndexFileName() const;
	string getSequenceFileName() const;
	string getHeaderFileName() const;
	int * getSequenceOffsetTable() const;
	int * getHeaderOffsetTable() const;
	int getNumberOfSequences() const;
	
	void indexInformation();
	void printSequenceName (int * places, int * scores);
};

#endif
