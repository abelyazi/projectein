#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "Database.h"

Database::Database(){
	index = "newE.fasta.pin";
	header = "./newE.fasta.phr";
	sequence = "./newE.fasta.psq";
	numberOfSequences = 0;
}


vector<uint32_t> Database::getSequenceOffsetTable(){
	return sequenceOffsetTable;
}

void Database::headerInformation(){

	uint32_t version;
	uint32_t dbType;
	uint32_t titleLength;
	string title;
	uint32_t timestampLength;
	string timestamp;
	//uint32_t numberOfSequences;
	uint64_t residueCount;
	uint32_t maximumSequence;
	
	ifstream pin;
	pin.open("newE.fasta.pin", ios::binary );
	pin.read((char*)&version, sizeof(uint32_t));
	version = __bswap_32(version);
	pin.read((char*)&dbType, sizeof(uint32_t));
	dbType= __bswap_32(dbType);
	//title
	pin.read((char*)&titleLength, sizeof(uint32_t));
	titleLength =__bswap_32(titleLength);
	char *titleSave = new char[titleLength];
	pin.read(titleSave, sizeof(char)*titleLength);
	titleSave[titleLength] = 0;
	title = titleSave;
	//timestamp
	pin.read((char*)&timestampLength, sizeof(uint32_t));
	timestampLength= __bswap_32(timestampLength);
	char *timestampSave = new char[timestampLength]; //a changer
	pin.read(timestampSave, sizeof(char)*timestampLength);
	timestampSave[timestampLength] = 0;
	timestamp = timestampSave;
	//number sequence
	pin.read((char*)&numberOfSequences, sizeof(uint32_t));
	numberOfSequences =__bswap_32(numberOfSequences);
	//residue count
	pin.read((char*)&residueCount, sizeof(uint64_t));
	//maximum sequence
	pin.read((char*)&maximumSequence, sizeof(uint32_t));
	maximumSequence =__bswap_32(maximumSequence);
	//offsets
	headerOffsetTable = vector<uint32_t>(numberOfSequences+1);
	sequenceOffsetTable = vector<uint32_t>(numberOfSequences+1);
	pin.read((char*)&headerOffsetTable[0], sizeof(uint32_t)*(numberOfSequences+1));
	pin.read((char*)&sequenceOffsetTable[0], sizeof(uint32_t)*(numberOfSequences+1));
	for(unsigned int i = 0; i < numberOfSequences+1; i++){
		headerOffsetTable[i] = __bswap_32(headerOffsetTable[i]);
		sequenceOffsetTable[i] = __bswap_32(sequenceOffsetTable[i]);
	}
	pin.close();
}



void Database::printSequenceName (int place){

	// calculer l'offset (début et fin) du header file grâce à l'index
	int offset1 = headerOffsetTable[place];
	int offset2 = headerOffsetTable[place+1];
	
	// chercher le nom dans le header file grâce à l'offset
	ifstream header;
	header.open("./newE.fasta.phr");
	header.seekg(offset1);
    int length = offset2-offset1;
    char* buffer = new char[length];
    header.read(buffer,length);
    header.close();
    
    // impression du nom trouvé dans le header
    cout << "Nom trouvé dans le header = ";
    cout.write (buffer,length) << endl;
    
    delete[] buffer; 
    
}









