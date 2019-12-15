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
public:
	Algorithm(int R, int Q);
	void calculate(Sequence * s, Database * d, int m_matrixBlosum[28][28]);
	int calculateScore(vector<int> sequence, int tailleSequence, int tailleSequence2, int m_matrixBlosum[28][28], int * databaseSequence);
};

#endif



