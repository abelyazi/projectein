#ifndef ALGORITHM_H
#define ALGORITHM_H

#include "Sequence.h"
#include "Database.h" 

using namespace std;

class Algorithm
{
private:

	int penalExtension; // R
	int penalOpen; // Q
	
	int (*m_matrixBlosum)[28];
	
public:

	Algorithm(int R, int Q, int (*mB)[28], ofstream * pfichier);
	void calculate(Sequence * s, Database * d);
	int calculateScore(int * sequence, int tailleSequence, int tailleSequence2, int * databaseSequence);
	
};

#endif



