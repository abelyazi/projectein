#include "Sequence.h"

Sequence::Sequence(){
}

Sequence::Sequence(string s){
	ifstream unknown_sequence(s);
	string real_sequence;
	if (unknown_sequence.is_open()){
		
		getline(unknown_sequence,name);
		string line = "";
		while(getline(unknown_sequence,line))
		{
		real_sequence += line;
		}
		unknown_sequence.close();
	}
	
	map<char, int> conversion = { {'-', 0}, {'A', 1}, {'B', 2}, {'C', 3}, {'D', 4}, {'E', 5}, {'F', 6}, {'G', 7}, {'H', 8}, {'I', 9}, {'J', 27}, {'K', 10}, {'L', 11}, {'M', 12}, {'N', 13}, {'O', 26}, {'P', 14}, {'Q', 15}, {'R', 16}, {'S', 17}, {'T', 18}, {'U', 24}, {'V', 19}, {'W', 20}, {'X', 21}, {'Y', 22}, {'Z', 23}, {'*', 25} };
	
	// convertion de la séquence inconnue en une liste d'entiers	
	size = real_sequence.size()+1;
	for ( int i=0;i<real_sequence.size()+1;i++){
		sequence.push_back(conversion[real_sequence[i]]);
		
	}
}

string Sequence::getName(){
	return name;
}

vector<int> Sequence::getSequence(){
	return sequence;
}

int Sequence::getSize(){
	return size;
}

void Sequence::printSequence(){
	cout << "Séquence de base = ";
	for (int n : sequence){
		cout << n;
	}
	cout << endl;
	cout << "Nom de base = " << name << endl;
}









