#include "Sequence.h"

Sequence::Sequence(string s, ofstream * pfichier){
	
	// récupérer les informations de la séquence dans le sequence file
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
	} else {
		cout << "You do not have the right query sequence file name. The name of the file should be " << s << endl;
	}
	
	map<char, int> conversion = { {'-', 0}, {'A', 1}, {'B', 2}, {'C', 3}, {'D', 4}, {'E', 5}, {'F', 6}, {'G', 7}, {'H', 8}, {'I', 9}, {'J', 27}, {'K', 10}, {'L', 11}, {'M', 12}, {'N', 13}, {'O', 26}, {'P', 14}, {'Q', 15}, {'R', 16}, {'S', 17}, {'T', 18}, {'U', 24}, {'V', 19}, {'W', 20}, {'X', 21}, {'Y', 22}, {'Z', 23}, {'*', 25} };
	
	// convertion de la séquence inconnue en une liste d'entiers	
	size = real_sequence.size()+1;
	sequence = new int[size]; //
	for ( int i=0;i<real_sequence.size()+1;i++){
		sequence[i] = (conversion[real_sequence[i]]);
	}
	
	// écrire les informations de la séquence de base dans le fichier texte
	if (pfichier->is_open()){
		*pfichier << "PROTEIN INFORMATION :" << endl;
		*pfichier << "Query sequence file : " << s << endl;
		*pfichier << "Query sequence size : " << size-1 << endl;
		*pfichier << "Query sequence name (for the verification) : " << name << endl << endl;
	}	
}

string Sequence::getName() const{
	return name;
}

int* Sequence::getSequence() const{
	return sequence;
}

int Sequence::getSize() const{
	return size;
}

void Sequence::printSequence(){
	cout << "Séquence de base = ";
	for (int i=0; i<size ; i++){
		cout << sequence[i];
	}
	cout << endl;
}









