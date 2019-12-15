#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <cstring>
#include <sstream>

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
	void calculate(Sequence * s, Database * d);
};

#endif




