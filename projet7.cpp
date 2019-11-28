#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include<iostream>
#include<fstream>
using namespace std;

int main(int argc, char **argv){
	
	// récupération de la sequence de base
	if (argc < 2){
		cout << "Aucune séquence inconnue fournie" << endl;
		return 1;
	}
	ifstream unknown_sequence(argv[1]); // P00533.fasta
	string real_name, real_sequence;
	if (unknown_sequence.is_open()){
		getline(unknown_sequence,real_name);
		getline(unknown_sequence,real_sequence);
		unknown_sequence.close();
	}
	
	// table de convertion psq
	vector<int> convertionTable = {1,2,3,4,5,6,7,8,9,27,10,11,12,13,26,14,15,16,17,18,24,19,20,21,22,23};
	
	// convertion de la séquence inconnue en une liste d'entiers	
	vector<int> sequence;
	for (int i=0;i<real_sequence.size();i++){
		sequence.push_back(convertionTable[real_sequence[i] - 65]);
	}
	
	// impression de la séquence de base et de son nom
	cout << "Séquence de base = ";
	for (int n : sequence){
		cout << n;
	}
	cout << endl;
	cout << "Nom de base = " << real_name << endl;

	// Chercher les informations nécessaires dans l'index
	
	uint32_t version;
	uint32_t dbType;
	uint32_t titleLength;
	string title;
	uint32_t timestampLength;
	string timestamp;
	uint32_t numberOfSequences;
	uint64_t residueCount;
	uint32_t maximumSequence;
	
	ifstream pin;
	pin.open("uniprot_sprot.fasta.pin", ios::binary );
	pin.read((char*)&version, sizeof(uint32_t));
	version = __bswap_32(version);
	pin.read((char*)&dbType, sizeof(uint32_t));
	dbType= __bswap_32(dbType);
	//title
	pin.read((char*)&titleLength, sizeof(uint32_t));
	titleLength =__bswap_32(titleLength);
	char titleSave[titleLength];
	pin.read((char*)titleSave, sizeof(char)*titleLength);
	titleSave[titleLength] = 0;
	title = titleSave;
	//timestamp
	pin.read((char*)&timestampLength, sizeof(uint32_t));
	timestampLength= __bswap_32(timestampLength);
	char timestampSave[timestampLength]; //a changer
	pin.read((char*)&timestampSave, sizeof(char)*timestampLength);
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
	//cout << version << endl << dbType << endl << titleLength << endl << title << timestampLength << endl <<timestamp << endl<< numberOfSequences << endl << residueCount << endl << maximumSequence << endl;
	//offsets
	vector<uint32_t> headerOffsetTable = vector<uint32_t>(numberOfSequences+1);
	vector<uint32_t> sequenceOffsetTable = vector<uint32_t>(numberOfSequences+1);
	pin.read((char*)&headerOffsetTable[0], sizeof(uint32_t)*(numberOfSequences+1));
	pin.read((char*)&sequenceOffsetTable[0], sizeof(uint32_t)*(numberOfSequences+1));
	for(int i = 0; i < numberOfSequences+1; i++){
		headerOffsetTable[i] = __bswap_32(headerOffsetTable[i]);
		sequenceOffsetTable[i] = __bswap_32(sequenceOffsetTable[i]);
	}
	pin.close();

	
	// recherche dans la database de la séquence
	ifstream database;
	database.open("./uniprot_sprot.fasta.psq");
	vector<int> foundSequence(sequence.size(),0);
	int place;
	for (int p = 0; p<sequenceOffsetTable.size(); p++){
		database.seekg(sequenceOffsetTable[p]);
		//if (sequence.size() == sequenceOffsetTable[p+1]-sequenceOffsetTable[p]){
			for (int i=0; i<sequence.size(); i++){
				char x;
				//database.seekg(sequenceOffsetTable[p]);
				if (database.get(x)){
					if (int(x) == sequence[i]){
						foundSequence[i] = int(x);
						if (i==sequence.size()-1){
							place = p;
							p = sequenceOffsetTable.size();
						}
					}
					else {
						i = sequence.size();
					}
				}
			}
		//}
	}
	
	// impression de la séquence trouvée dans la database
	cout << "Séquence trouvée dans la database = ";
	for (int n:foundSequence){
		cout << n;
	} cout << endl;
	
	// chercher les offset dans le header file grâce à l'index
	int offset1 = headerOffsetTable[place];
	int offset2 = headerOffsetTable[place+1];
	
	// chercher le nom dans le header file grâce à l'offset
	ifstream header;
	header.open("./uniprot_sprot.fasta.phr");
	header.seekg(offset1);
    int length = offset2-offset1;
    char* buffer = new char[length];
    header.read(buffer,length);
    header.close();
    
    // impression du nom trouvé dans le header
    cout << "Nom trouvé dans le header = ";
    cout.write (buffer,length) << endl;

	return 0;
	
}

