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

int main(int argc, char **argv)
{
	Database * db = new Database();
	if (argc < 2){
		cout << "Aucune séquence inconnue fournie" << endl;
		return 1;
	}
	Sequence * seq = new Sequence(argv[1]);
	Algorithm * algo = new Algorithm(1, 12);
	algo->calculate(seq, db);
	
	delete db;
	delete seq;
	delete algo;

	return 0;	
}


