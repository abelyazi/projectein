#include "Database.h"

Database::Database(string n, ofstream * pf){
	name = n;
	pfichier = pf;
	indexInformation();
}

string Database::getName() const{
	return name;
}

vector<uint32_t> Database::getSequenceOffsetTable() const{
	return sequenceOffsetTable;
}

vector<uint32_t> Database::getHeaderOffsetTable() const{
	return headerOffsetTable;
}

int Database::getNumberOfSequences() const{
	return numberOfSequences;
}

void Database::indexInformation(){ // récupérer les informations dans l'index

	uint32_t version;
	uint32_t dbType;
	uint32_t titleLength;
	string title;
	uint32_t timestampLength;
	string timestamp;
	uint64_t residueCount;
	uint32_t maximumSequence;
	
	ifstream indexFile;
	indexFile.open(name+".pin", ios::binary );
	if (indexFile.is_open()){
		indexFile.read((char*)&version, sizeof(uint32_t));
		version = __bswap_32(version);
		indexFile.read((char*)&dbType, sizeof(uint32_t));
		dbType= __bswap_32(dbType);
		//title
		indexFile.read((char*)&titleLength, sizeof(uint32_t));
		titleLength =__bswap_32(titleLength);
		char *titleSave = new char[titleLength];
		indexFile.read(titleSave, sizeof(char)*titleLength);
		titleSave[titleLength] = 0;
		title = titleSave;
		//timestamp
		indexFile.read((char*)&timestampLength, sizeof(uint32_t));
		timestampLength= __bswap_32(timestampLength);
		char *timestampSave = new char[timestampLength]; //a changer
		indexFile.read(timestampSave, sizeof(char)*timestampLength);
		timestampSave[timestampLength] = 0;
		timestamp = timestampSave;
		//number sequence
		indexFile.read((char*)&numberOfSequences, sizeof(uint32_t));
		numberOfSequences =__bswap_32(numberOfSequences);
		//residue count
		indexFile.read((char*)&residueCount, sizeof(uint64_t));
		//maximum sequence
		indexFile.read((char*)&maximumSequence, sizeof(uint32_t));
		maximumSequence =__bswap_32(maximumSequence);
		//offsets
		headerOffsetTable = vector<uint32_t>(numberOfSequences+1);
		sequenceOffsetTable = vector<uint32_t>(numberOfSequences+1);
		indexFile.read((char*)&headerOffsetTable[0], sizeof(uint32_t)*(numberOfSequences+1));
		indexFile.read((char*)&sequenceOffsetTable[0], sizeof(uint32_t)*(numberOfSequences+1));
		for(unsigned int i = 0; i < numberOfSequences+1; i++){
			headerOffsetTable[i] = __bswap_32(headerOffsetTable[i]);
			sequenceOffsetTable[i] = __bswap_32(sequenceOffsetTable[i]);
		}
		indexFile.close();
	} else {
		cout << "You do not have the right index file name. The name of the file should be " << name+".pin" << endl;
	}
	
	
	// écrire les informations de la database dans le fichier texte
	if (pfichier->is_open()){
		*pfichier << "DATABASE INFORMATION :" << endl;
		*pfichier << "Version : " << version << endl;
		*pfichier << "Type : " << dbType << endl;
		*pfichier << "Title : " << title << endl;
		*pfichier << "Database Time : " << timestamp << endl;
		*pfichier << "Number of sequences : " << numberOfSequences << endl;
		*pfichier << "Number of residues : " << residueCount << endl;
		*pfichier << "Longest db seq : " << maximumSequence << " residues" << endl << endl;
	}

}


void Database::printSequenceName (int * places, int * scores){
	
	*pfichier << "DISPLAY OF THE 50 BEST ALIGNEMENTS : " << endl << endl;
	
	// ouvrir le header file
	ifstream headerFile;
	headerFile.open(name+".phr");
	
	if (headerFile.is_open()){
		for (int i = 49 ; i>=0; i-- ){
			*pfichier << "gnl|BL_ORD_ID|" << places[i] << " ";
			
			// trouver les offsets de la séquence
			int offset1 = headerOffsetTable[places[i]];
			int offset2 = headerOffsetTable[places[i]+1];
			
			// chercher le nom dans le header file grâce à l'offset
			headerFile.seekg(offset1);
			int length = offset2-offset1+10;
			char* buffer = new char[length];
			headerFile.read(buffer,length);
			
			// impression du nom trouvé dans le header
			bool read = false;
			for (int j=0; j<53 ;j++){
				if ((unsigned int)(unsigned char)(buffer[j]) == 0x73){
					read = true;
				}
				if (read){
					*pfichier << buffer[j];
				}
			} *pfichier << "...  " << int((0.267*scores[i]+3.34)/0.6931471806) << endl;
			delete[] buffer; 
		}
	} else {
		cout << "You do not have the right header file name. The name of the file should be " << name+".phr" << endl;
	}

	headerFile.close();
}
