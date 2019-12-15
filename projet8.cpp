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

#include "Sequence.h"
#include "Database.h"
#include "Algorithm.h"

using namespace std;
void blosumCreation(int m_matrixBlosum[28][28]);

int main(int argc, char **argv)
{
	int m_matrixBlosum[28][28];
	blosumCreation(m_matrixBlosum);
	
	Database * database = new Database();
	if (argc < 2){
		cout << "Aucune séquence inconnue fournie" << endl;
		return 1;
	}
	
	Sequence * sequence = new Sequence(argv[1]);
	
	Algorithm * algorithm = new Algorithm(1, 12);
	algorithm->calculate(sequence, database, m_matrixBlosum);
	
	delete database;
	delete sequence;
	delete algorithm;

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


