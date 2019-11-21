#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include<iostream>
#include<fstream>
using namespace std;

int main(){
	
	// récupération de la sequence de base
	ifstream unknown_sequence("./P00533.fasta");
	string real_name, real_sequence;
	if (unknown_sequence.is_open()){
		getline(unknown_sequence,real_name);
		getline(unknown_sequence,real_sequence);
		unknown_sequence.close();
	}
	
	// table de convertion psq
	vector<int> convertion_table = {1,2,3,4,5,6,7,8,9,27,10,11,12,13,26,14,15,16,17,18,24,19,20,21,22,23};
	
	// convertion de la séquence inconnue en une liste d'entiers	
	vector<int> s;
	for (int i=0;i<real_sequence.size();i++){
		int n = real_sequence[i] - 65;
		s.push_back(convertion_table[n]);
	}
	
	// impression de la séquence de base
	cout << "séquence de base = " << endl;
	for (int a : s){
		cout << a << " ";
	}
	cout << endl;
	
	// recherche dans la database
	ifstream database("./uniprot_sprot.fasta.psq", ios::binary); // , ios::bin
	
	char x;
	int i = 0;
	vector<int> db_s(s.size(),0);
	int y;

	if(database.is_open()) {
		while(database.get(x)){ 
			y = int(x);
			if (y==s[i]){
				db_s[i] = y;
				i++;
				if (i>=s.size()){
					break;
				}
			}
			else {
				i = 0;
			}
		}
	}
	
	// impression de la séquence trouvée dans la database
	cout << "séquence trouvée dans la database = " << endl;
	for (int b:db_s){
		cout << b << " ";
	}
	cout << endl;
	
	
	return 0;
	
}
