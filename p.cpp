/* 
Ammi Haroun
Belyazid Ali
de Wouters Louise
Remise intermédiaire - 29/11/2019
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
 //typedef int Array2x2[28][28];
struct resultat{
	int score;
	int id;
} ;

struct thread_args{
    int depart;
    int fin;
    vector<uint32_t>* sequenceOffsetTable;
    vector<int>* sequence;
    int* databaseTable;
    int tailleSequence;
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
	
	vector<int> convertionTable = {1,2,3,4,5,6,7,8,9,27,10,11,12,13,26,14,15,16,17,18,24,19,20,21,22,23};
	map<char, int> conversion = { {'-', 0}, {'A', 1}, {'B', 2}, {'C', 3}, {'D', 4}, {'E', 5}, {'F', 6}, {'G', 7}, {'H', 8}, {'I', 9}, {'J', 27}, {'K', 10}, {'L', 11}, {'M', 12}, {'N', 13}, {'O', 26}, {'P', 14}, {'Q', 15}, {'R', 16}, {'S', 17}, {'T', 18}, {'U', 24}, {'V', 19}, {'W', 20}, {'X', 21}, {'Y', 22}, {'Z', 23}, {'*', 25} };
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
		cout << "Aucune séquence inconnue fournie" << endl;
		return 1;
	}
	ifstream unknown_sequence(argv[1]); // P00533.fasta
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
	
	// table de convertion psq
	//vector<int> convertionTable = {1,2,3,4,5,6,7,8,9,27,10,11,12,13,26,14,15,16,17,18,24,19,20,21,22,23};
	
	// convertion de la séquence inconnue en une liste d'entiers	
	vector<int> sequence;
	int tailleSequence = real_sequence.size()+1;
	
	for ( int i=0;i<real_sequence.size()+1;i++){
		sequence.push_back(conversion[real_sequence[i]]);
		
	}
	
	// impression de la séquence de base et de son nom (pour la vérification)
	/*cout << "Séquence de base = ";
	for (int n : sequence){
		cout << n;
	}
	cout << endl;
	cout << "Nom de base = " << real_name << endl << endl; */

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
	//cout << version << endl << dbType << endl << titleLength << endl << title << timestampLength << endl <<timestamp << endl<< numberOfSequences << endl << residueCount << endl << maximumSequence << endl;
	//offsets
	vector<uint32_t> headerOffsetTable = vector<uint32_t>(numberOfSequences+1);
	vector<uint32_t> sequenceOffsetTable = vector<uint32_t>(numberOfSequences+1);
	pin.read((char*)&headerOffsetTable[0], sizeof(uint32_t)*(numberOfSequences+1));
	pin.read((char*)&sequenceOffsetTable[0], sizeof(uint32_t)*(numberOfSequences+1));
	for(unsigned int i = 0; i < numberOfSequences+1; i++){
		headerOffsetTable[i] = __bswap_32(headerOffsetTable[i]);
		sequenceOffsetTable[i] = __bswap_32(sequenceOffsetTable[i]);
	}
	pin.close();

	// recherche de la séquence dans le sequence file
	pthread_t new_thread;
	pthread_t new_thread2;
	ifstream database;
	database.open("./newE.fasta.psq");
	int* databaseTable = new int[sequenceOffsetTable[numberOfSequences]];
	//cout<<"XDDDDDDDd"<<endl;
	
	for (unsigned int i=0; i<sequenceOffsetTable[numberOfSequences]; i++){
		char x;
		if (database.get(x)){
			databaseTable[i] = (x);
			//cout<< databaseTable[i];
		}
	}
	resultat* result= new resultat[numberOfSequences];
	/*int* score = new int[numberOfSequences];
	int* id = new int[numberOfSequences];
	result->score = score;
	result->id = id;*/
	//int a = sequence.size();
	//cout <<a<<endl;
	//cout<<"OHLALALALA"<<endl;
	//vector<int> foundSequence(sequence.size());
	int place = 0;
	int penalExtension = 1;
	int penalOpen = 11 + penalExtension;
	
	int bp = 3;
	int scoreMaxdb = 0;
	//cout<<"azlezla"<<endl;
	
	thread_args *args = new thread_args;
	thread_args *args2 = new thread_args;
		int depart = 0;
		int fin = 7000;
		int depart2 = 3500;
		int fin2 = 7000;
		int scoreMaxdb2 = 0;
		if (args != NULL)
		{
			args->depart = depart;
			
			args->fin = fin;
			args->sequenceOffsetTable = &sequenceOffsetTable;
			args->databaseTable = databaseTable;
			args->tailleSequence = tailleSequence;
			args->sequence = &sequence;
			args->penalExtension = &penalExtension;
			args->penalOpen = &penalOpen;
			args->m_matrixBlosum = &m_matrixBlosum;
			args->scoreMaxdb = &scoreMaxdb;
			args->resultats = result;
				//cout<<tailleSequence<<" Parfait"<<endl;

			pthread_create (&new_thread, NULL, maxscoreId, (void *)args);
		}
		//cout<<"lallalll "<<args->depart<<endl;
		/*		if (args2 != NULL)
		{
			args2->depart = depart2;
			args2->fin = fin2;
			args2->sequenceOffsetTable = &sequenceOffsetTable;
			args2->databaseTable = databaseTable;
			args2->tailleSequence = tailleSequence;
			args2->sequence = &sequence;
			args2->penalExtension = &penalExtension;
			args2->penalOpen = &penalOpen;
			args2->m_matrixBlosum = &m_matrixBlosum;
			args2->scoreMaxdb = &scoreMaxdb2;
								//cout<<tailleSequence<<" Parfait 2"<<endl;
			pthread_create (&new_thread2, NULL, maxscoreId, (void *)args2);
		}
		pthread_join (new_thread2, NULL);*/
		pthread_join (new_thread, NULL);
		//cout<<"non "<<scoreMaxdb<<" Ayaa"<<endl;
		
		/*if(scoreMaxdb2 > scoreMaxdb)
		{scoreMaxdb = scoreMaxdb2;}*/
		
		sort(result,result+numberOfSequences,maxStruct);
		cout<<" "<<result[0].score<<" "<<result[1].score<<endl; 

		
	//cout << "Best score db : " << scoreMaxdb <<endl;
	database.close();
	/*
	// impression de la séquence trouvée dans le sequence file
	cout << "Séquence trouvée dans la database = ";
	for (int n:foundSequence){
		cout << n;
	} cout << endl;*/
	
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
    cout << "best score :" << scoreMaxdb << " place :" <<place<< endl;
    cout << "Nom trouvé dans le header = ";
    cout.write (buffer,length) << endl;
    
    delete[] buffer; 
    delete[] databaseTable;
    delete args;

	return 0;	
}

void *maxscoreId (void *data)
{
	//int scoreMaxdb = 0;
	thread_args *args = (thread_args*)data;
	//printf("nombre 1 : %d nombre 2 : %d    \n",*(args->penalOpen),*(args->penalExtension));
	

	for (unsigned int p = args->depart; p<args->sequenceOffsetTable->size()-1 && p<args->fin; p++){
		//cout<<"NOOOOOOOOOON"<<endl;
		int tailleSequence2 = (*(args->sequenceOffsetTable))[p+1]-(*(args->sequenceOffsetTable))[p]-1;
		
		//cout<<"ICCCIIIIIIIii"<<endl;
		int* shady = new int[tailleSequence2];
		//cout<<"hassan est moche"<<endl;
		for (unsigned int i=0; i<tailleSequence2; i++){
		
				shady[i] = args->databaseTable[((*(args->sequenceOffsetTable))[p])+i];
				//cout<< args->databaseTable[((*(args->sequenceOffsetTable))[p])+i];
			
		}
		//cout<<"fin "<<endl;
		int* matrixH = new int[args->tailleSequence+1];
		int* matrixE = new int[args->tailleSequence+1];
		int F = 0;
		int Hhaut = 0;
		int Hdiag = 0;
		int Hmaintenant = 0;
		int scoreMaxseq =0;
		
		for (unsigned int i=0; i< args->tailleSequence+1;i++)
			{
			matrixE[i] = 0;
			matrixH[i] = 0;
			
			}
					
		//cout << "taille i :" << tailleSequence << " taille j :" << tailleSequence2 << endl;
		for (unsigned int j=1; j<tailleSequence2+1; j++) { // ligne
			for (unsigned int i=1; i<args->tailleSequence+1; i++) { //colonne
			//Eij
				// cout<<i<<"  "<<j<<endl;
				matrixE[i]=max(matrixH[i]-*(args->penalOpen),matrixE[i]-*(args->penalExtension)); //  ça c'est logique
				F=max(Hhaut-*(args->penalOpen), F-*(args->penalExtension));
				Hdiag = matrixH[i-1] + (*(args->m_matrixBlosum))[shady[j-1]][(*(args->sequence))[i-1]];
				Hmaintenant = max(Hdiag, matrixE[i]);
				Hmaintenant = max(Hmaintenant, F);
				Hmaintenant = max(Hmaintenant, 0);

				if(Hmaintenant > scoreMaxseq){
					scoreMaxseq = Hmaintenant;
				}
				matrixH[i-1] = Hhaut;
				//Eleft[i] = E;
				Hhaut= Hmaintenant;
			}
		
			matrixH[args->tailleSequence] = Hmaintenant;
			
			Hhaut = 0;//printf("\n");
		}
		//cout<< p << "   " << scoreMaxseq<<endl;
		//printf("%d %d\n",matrixH[coord[0]*2 +coord[1]],p);

		//cout<< scoreMaxseq<<endl;
		if(scoreMaxseq > *(args->scoreMaxdb))	*(args->scoreMaxdb) = scoreMaxseq;
		(args->resultats[p]).score = scoreMaxseq;
		(args->resultats[p]).id = p;
		delete[] matrixH;
		delete[] matrixE;
		delete[] shady;
	}
	

}
