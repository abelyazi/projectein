/* 
Ammi Haroun
Belyazid Ali
de Wouters Louise
code annexe : utilisation de thread pour optimiser
*/

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <cstring>
#include <sstream>
#include <pthread.h>
#include <sys/sysinfo.h>
#include <algorithm>
#define max(x,y) ( x < y ? y : x ) 

using namespace std;
struct resultat{
	int score;
	int id;
} ;

struct thread_args{
    int start;
    int end;
    vector<uint32_t>* sequenceOffsetTable;
    vector<int>* querySequence;
    int* databaseTable;
    int querySequenceSize;
    int* penalExtension;
    int* penalOpen;
    int (*m_matrixBlosum)[28][28];
    int* scoreMaxdb; 
    resultat* resultats;
    
}  ;
void * maxscoreId (void *data);
bool maxStruct (resultat i,resultat j) { return (i.score>j.score); }
int main(int argc, char **argv)
{

	// table de convertion
	map<char, int> conversion = { {'-', 0}, {'A', 1}, {'B', 2}, {'C', 3}, {'D', 4}, {'E', 5}, {'F', 6}, {'G', 7}, {'H', 8}, {'I', 9}, {'J', 27}, {'K', 10}, {'L', 11}, {'M', 12}, {'N', 13}, {'O', 26}, {'P', 14}, {'Q', 15}, {'R', 16}, {'S', 17}, {'T', 18}, {'U', 24}, {'V', 19}, {'W', 20}, {'X', 21}, {'Y', 22}, {'Z', 23}, {'*', 25} };
	
	// création de la matrice BLOSUM
	ifstream m_blosumName("BLOSUM62");
	int m_matrixBlosum[28][28];
	for (int i = 0; i<28; i++) {
		for (int j = 0; j<28; j++) {
			m_matrixBlosum[i][j] = 0; 
		}
	}
	if (m_blosumName)  
	{
		string s = "";
		map<int, int> conversionChar;
		int j = 1;
		while (getline( m_blosumName, s )) 
		{	
			if (s[0] != '#')
			{		
				if(s[0]==' ')
				{
					for (int i = 0; i < s.size(); i++) 
					{
						if (s[i] != ' ') 
						{			
							conversionChar.insert({i/3, conversion[s[i]]});
						}
					}
				}
				else
				{
					for (int i = 3; i < 73; i+=3)
					{
						if(s[i-1] == '1') 
						m_matrixBlosum[conversionChar[i/3]][conversionChar[j]]= (10 + (int) s[i]-48) *((s[i-1]=='-' ) ? -1 : 1) ;
						else m_matrixBlosum[conversionChar[i/3]][conversionChar[j]]= ((int) s[i]-48) *((s[i-1]=='-' ) ? -1 : 1) ;
					} 
					j++;
				}
			}
		}

	}
	
	// récupération de la sequence de base
	if (argc < 2){
		cout << "Argument missing" << endl;
		return 1;
	}
	ifstream unknown_sequence(argv[1]);
	string real_name, real_sequence;
	if (unknown_sequence.is_open()){
		getline(unknown_sequence,real_name);
		string line = "";
		while(getline(unknown_sequence,line))
		{
		real_sequence += line;
		}
		unknown_sequence.close();
	}
	
	// convertion de la séquence inconnue en une liste d'entiers	
	vector<int> querySequence;
	int querySequenceSize = real_sequence.size()+1;
	for ( int i=0;i<real_sequence.size()+1;i++){
		querySequence.push_back(conversion[real_sequence[i]]);
	}

	// chercher les informations nécessaires dans l'index
	uint32_t version;
	uint32_t dbType;
	uint32_t titleLength;
	string title;
	uint32_t timestampLength;
	string timestamp;
	uint32_t numberOfSequences;
	uint64_t residueCount;
	uint32_t maximumSequence;
	
	ifstream indexFile;
	indexFile.open("database.fasta.pin", ios::binary );
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
	vector<uint32_t> headerOffsetTable = vector<uint32_t>(numberOfSequences+1);
	vector<uint32_t> sequenceOffsetTable = vector<uint32_t>(numberOfSequences+1);
	indexFile.read((char*)&headerOffsetTable[0], sizeof(uint32_t)*(numberOfSequences+1));
	indexFile.read((char*)&sequenceOffsetTable[0], sizeof(uint32_t)*(numberOfSequences+1));
	for(unsigned int i = 0; i < numberOfSequences+1; i++){
		headerOffsetTable[i] = __bswap_32(headerOffsetTable[i]);
		sequenceOffsetTable[i] = __bswap_32(sequenceOffsetTable[i]);
	}
	indexFile.close();

	// recherche de la séquence dans le sequence file
	
	// ouvrir le sequence file
	ifstream database;
	database.open("./database.fasta.psq");
	
	// mettre la database dans un tableau
	int* databaseTable = new int[sequenceOffsetTable[numberOfSequences]];
	for (unsigned int i=0; i<sequenceOffsetTable[numberOfSequences]; i++){
		char x;
		if (database.get(x)){
			databaseTable[i] = (x);
		}
	}
	
	int place = 0;
	int penalExtension = 1;
	int penalOpen = 11 + penalExtension;
	//int bp = 3;
	int scoreMaxdb = 0;
	
	resultat* result= new resultat[numberOfSequences];
	
	// créer les thread
	pthread_t new_thread;
	pthread_t new_thread2;
	thread_args *args = new thread_args;
	thread_args *args2 = new thread_args;
	int start1 = 0;
	int end1 = 3500;
	int start2 = 3501;
	int end2 = 7000;
	int scoreMaxdb2 = 0;
	if (args != NULL)
	{
		args->start = start1;
		args->end = end1;
		args->sequenceOffsetTable = &sequenceOffsetTable;
		args->databaseTable = databaseTable;
		args->querySequenceSize = querySequenceSize;
		args->querySequence = &querySequence;
		args->penalExtension = &penalExtension;
		args->penalOpen = &penalOpen;
		args->m_matrixBlosum = &m_matrixBlosum;
		args->scoreMaxdb = &scoreMaxdb;
		args->resultats = result;

		pthread_create (&new_thread, NULL, maxscoreId, (void *)args);
	}	
	if (args2 != NULL)
	{
		args2->start = start2;
		args2->end = end2;
		args2->sequenceOffsetTable = &sequenceOffsetTable;
		args2->databaseTable = databaseTable;
		args2->querySequenceSize = querySequenceSize;
		args2->querySequence = &querySequence;
		args2->penalExtension = &penalExtension;
		args2->penalOpen = &penalOpen;
		args2->m_matrixBlosum = &m_matrixBlosum;
		args2->scoreMaxdb = &scoreMaxdb2;
		args2->resultats = result;
		
		pthread_create (&new_thread2, NULL, maxscoreId, (void *)args2);
	}
	
	pthread_join (new_thread2, NULL);
	pthread_join (new_thread, NULL);
	
	sort(result,result+numberOfSequences,maxStruct);
	
	// affichage des 10 meilleurs scores
	cout << endl << "Display of the 10 best alignements" << endl << endl;
	for (int i = 0; i<10 ;i++){
		cout << "score " << i << " : " << result[i].score << " à la place " << result[i].id << endl;
	}

	database.close();

	return 0;	
}

void *maxscoreId (void *data)
{
	thread_args *args = (thread_args*)data;
	
	for (unsigned int p = args->start; p<args->sequenceOffsetTable->size()-1 && p<args->end; p++){
		int databaseSequenceSize = (*(args->sequenceOffsetTable))[p+1]-(*(args->sequenceOffsetTable))[p]-1;
		
		int* databaseSequence = new int[databaseSequenceSize];
		
		for (unsigned int i=0; i<databaseSequenceSize; i++){
		
				databaseSequence[i] = args->databaseTable[((*(args->sequenceOffsetTable))[p])+i];
				
		}
		int* matrixH = new int[args->querySequenceSize+1];
		int* matrixE = new int[args->querySequenceSize+1];
		int F = 0;
		int H_up = 0;
		int H_diag = 0;
		int H = 0;
		int seqScoreMax =0;
		
		for (unsigned int i=0; i< args->querySequenceSize+1;i++)
			{
			matrixE[i] = 0;
			matrixH[i] = 0;
			
			}
		for (unsigned int j=1; j<databaseSequenceSize+1; j++) {
			for (unsigned int i=1; i<args->querySequenceSize+1; i++) {
				matrixE[i]=max(matrixH[i]-*(args->penalOpen),matrixE[i]-*(args->penalExtension));
				F=max(H_up-*(args->penalOpen), F-*(args->penalExtension));
				H_diag = matrixH[i-1] + (*(args->m_matrixBlosum))[databaseSequence[j-1]][(*(args->querySequence))[i-1]];
				H = max(H_diag, matrixE[i]);
				H = max(H, F);
				H = max(H, 0);

				if(H > seqScoreMax){
					seqScoreMax = H;
				}
				matrixH[i-1] = H_up;
				H_up= H;
			}
		
			matrixH[args->querySequenceSize] = H;
			
			H_up = 0;
		}

		if(seqScoreMax > *(args->scoreMaxdb))	*(args->scoreMaxdb) = seqScoreMax;
		(args->resultats[p]).score = seqScoreMax;
		(args->resultats[p]).id = p;
		delete[] matrixH;
		delete[] matrixE;
		delete[] databaseSequence;
	}
	

}
