/* 
Ammi Haroun
Belyazid Ali
de Wouters Louise
Projet cpa - 18/12/2019
*/

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <cstring>
#include <sstream>

#include "Sequence.h"
#include "Database.h"
#include "Algorithm.h"

using namespace std;
void blosumCreation(int m_matrixBlosum[28][28]);

int main()
{
	// générer la matrice blossum
	int m_matrixBlosum[28][28];
	blosumCreation(m_matrixBlosum);
	
	// ouvrir le fichier pour écrire les informations et les résultats
	ofstream fichier;
	fichier.open("./results.txt");
	if (fichier.is_open()){
		fichier << "Ammi Haroun - informatique" << endl << "Belyazid Ali - informatique" << endl << "de Wouters Louise - biomedical" << endl;
		fichier << "Projet cpa - remise le 18 décembre 2019" << endl << endl;
	}
	
	ofstream * pfichier = &fichier;
	
	// créer les trois objets utiles au programme
	Database * database = new Database(pfichier);

	Sequence * sequence = new Sequence("./proteine.fasta", pfichier);
	
	Algorithm * algorithm = new Algorithm(1, 12, m_matrixBlosum, pfichier);
	
	// lancer l'algorithme
	algorithm->calculate(sequence, database);
	
	delete database;
	delete sequence;
	delete algorithm;
	fichier.close();

	return 0;	
}

void blosumCreation(int m_matrixBlosum[28][28]){
	
	map<char, int> conversion = { {'-', 0}, {'A', 1}, {'B', 2}, {'C', 3}, {'D', 4}, {'E', 5}, {'F', 6}, {'G', 7}, {'H', 8}, {'I', 9}, {'J', 27}, {'K', 10}, {'L', 11}, {'M', 12}, {'N', 13}, {'O', 26}, {'P', 14}, {'Q', 15}, {'R', 16}, {'S', 17}, {'T', 18}, {'U', 24}, {'V', 19}, {'W', 20}, {'X', 21}, {'Y', 22}, {'Z', 23}, {'*', 25} };
	ifstream m_blosumName("BLOSUM62");
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
}


