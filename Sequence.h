#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <fstream>
#include <iostream>
#include <vector>
#include <map>

using namespace std;

class Sequence
{
private:

	string name;
	int * sequence;
	int size;
	
public:

	Sequence(string s, ofstream * pfichier);
	string getName() const;
	int * getSequence() const;
	int getSize() const;
	void printSequence();
};

#endif



