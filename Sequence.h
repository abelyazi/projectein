#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <cstring>
#include <sstream>
using namespace std;

class Sequence
{
private:
	string name;
	vector<int> sequence;
	int size;
public:
	Sequence();
	Sequence(string s);
	string getName();
	vector<int> getSequence();
	int getSize();
	void printSequence();
};

#endif



