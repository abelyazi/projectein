#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "Algorithm.h"

Algorithm::Algorithm(int R, int Q){
	penalExtension = R; // R
	penalOpen = Q; // Q
}

void Algorithm::calculate(Sequence * seq, Database * db){
	//MATRICE BLOSUM
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
	
	ifstream database;
	database.open("./newE.fasta.psq");
	//int penalExtension = 1; // R
	//int penalOpen = 11 + penalExtension; // Q
	vector<int> sequence = seq->getSequence(); //
	int tailleSequence = seq->getSize(); //
	db->headerInformation(); //
	vector<uint32_t> sequenceOffsetTable = db->getSequenceOffsetTable(); //
	vector<uint32_t> headerOffsetTable = db->getHeaderOffsetTable(); //
	
	int results[5] = {0};
	int scores[5] = {0};

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
		//int score = calculateScore(tailleSequence, tailleSequence2, databaseSequence, m_matrixBlosum);
		
		int* matrixH = new int[tailleSequence*tailleSequence2];
		int* matrixE = new int[tailleSequence*tailleSequence2];
		int* matrixF = new int[tailleSequence*tailleSequence2];
		int score = 0;
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
				
				//score
				score = matrixH[i*tailleSequence2 + j];
			}	
		}
		
		delete matrixH;
		delete matrixE;
		delete matrixF;
		
		
		if (score > scores[0]){
			scores[0] = score;
			results[0] = p;
			int i = 0;
			while (scores[i] > scores[i+1]){
				int temp = scores[i];
				scores[i] = scores[i+1];
				scores[i+1] = temp;
				int temp2 = results[i];
				results[i] = results[i+1];
				results[i+1] = temp2;
				i++;
				if (i == 4){
					break;
				}
			}
		}
	}
	
	database.close();
	int a = 0;
	for (int pl : results){
		cout << "score = " << scores[a] << endl;
		db->printSequenceName(pl);
		a++;
	}
	
}
