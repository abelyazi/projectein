/* 
Ammi Haroun
Belyazid Ali
de Wouters Louise
Projet - 18/12/2019
*/

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <cstring>
#include <sstream>

#include "Database.h"

using namespace std;

int main(int argc, char **argv)
{
	// générer la database 
	Database * db = new Database();
	
	// matrice BLOSUM
	
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
	
	//table de convertion psq
	vector<int> convertionTable = {1,2,3,4,5,6,7,8,9,27,10,11,12,13,26,14,15,16,17,18,24,19,20,21,22,23};
	
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

	// recherche de la séquence dans le sequence file
	ifstream database;
	database.open("./newE.fasta.psq");
	vector<int> foundSequence(sequence.size(),0);
	int place = 0;
	int penalExtension = 1;
	int penalOpen = 11 + penalExtension;
	int coord[2] = {0,0};
	int scoreLimite = 5000;
	
	vector<int> results;
	vector<int> scores;
	
	db->headerInformation();
	vector<uint32_t> sequenceOffsetTable = db->getSequenceOffsetTable();

	for (unsigned int p = 0; p<sequenceOffsetTable.size()-1; p++){
		database.seekg(sequenceOffsetTable[p]);
		int tailleSequence2 = sequenceOffsetTable[p+1]-sequenceOffsetTable[p];
		int databaseSequence[tailleSequence2] = {0};
		for (unsigned int i=0; i<tailleSequence2; i++){
			char x;
			if (database.get(x)){
				databaseSequence[i] = (x);
			}
		}
		int* matrixH = new int[tailleSequence*tailleSequence2];
		int* matrixE = new int[tailleSequence*tailleSequence2];
		int* matrixF = new int[tailleSequence*tailleSequence2];
		for (unsigned int i=0; i< tailleSequence*tailleSequence2;i++) {
			matrixE[i] = 0;
			matrixF[i] = 0;
			matrixH[i] = 0;
			}		
		for (unsigned int i=1; i<tailleSequence; i++) {
			for (unsigned int j=1; j<tailleSequence2; j++) {
				//Eij
				if(matrixE[i*tailleSequence2 + j-1] - penalExtension < matrixH[i*tailleSequence2 + j-1] - penalOpen ) {		
					matrixE[i*tailleSequence2 + j] = matrixH[i*tailleSequence2 + j-1] - penalOpen;
				}
				else {
					matrixE[i*tailleSequence2 + j] = matrixE[i*tailleSequence2 + j-1] - penalExtension;
				} 
				//FIJ 
				if(matrixF[(i-1)*tailleSequence2 + j] - penalExtension < matrixH[(i-1)*tailleSequence2 + j] - penalOpen ) {		
					matrixF[i*tailleSequence2 + j] = matrixH[(i-1)*tailleSequence2 + j] - penalOpen;
				}
				else {
					matrixF[i*tailleSequence2 + j] = matrixF[(i-1)*tailleSequence2 + j] - penalExtension;
				}
				//Hij
				matrixH[i*tailleSequence2 + j] = matrixH[(i-1)*tailleSequence2+j-1] + m_matrixBlosum[sequence[i-1]][databaseSequence[j-1]];
				if(matrixH[i*tailleSequence2 + j] < (matrixE[i*tailleSequence2 + j])) {
					matrixH[i*tailleSequence2 + j] = (matrixE[i*tailleSequence2 + j]);
				}
				if(matrixH[i*tailleSequence2 + j] < (matrixF[i*tailleSequence2 + j])) {
					matrixH[i*tailleSequence2 + j] = (matrixF[i*tailleSequence2 + j]);
				}
				if(matrixH[i*tailleSequence2 + j] < 0) {
					matrixH[i*tailleSequence2 + j] = 0;
				}
				//Scores
				if(matrixH[i*tailleSequence2 + j] > matrixH[coord[0]*tailleSequence2 +coord[1]]) {
					coord[0] = i;
					coord[1] = j;

				}
			}	
		}
		int score = matrixH[coord[0]*tailleSequence2 +coord[1]];
		if (score >= scoreLimite){
			results.push_back(p);
			scores.push_back(score);
			cout << "bon score de " << score << " à la place " << p << endl;
		}
		
		delete matrixH;
		delete matrixE;
		delete matrixF;
	}
	
	database.close();

	for (int pl : results){
		db->printSequenceName(pl);
	}

	delete db;

	return 0;	
}


