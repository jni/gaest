	/*
File:	test.cpp
Title:	Useful program to test the input of gaest, the est clustering program
		based on GA's
Descr:	invoke the program with the file containing the sequences as a command-
		line argument. The program will not check that the file exists
		or that the call was correct.
	*/


#include <iostream>
#include <vector>
#include <fstream>
#include "dna.h"
#include "dynamic.h"

int main( int argc, char** argv )
{
	vector<dna> sequences;
	dna temp;
	ifstream datafile( argv[1], ios::in );

	cerr	<< "File open. Reading in sequences" << endl;

	while( datafile >> temp )
	{
		cerr	<< "Sequence read." << endl;
		sequences.push_back(temp);
	}

	cerr	<< "Sequences are ready." << endl;
	cout	<< "Enter command: 1-print, 2-align." << endl;
	int c, i, j;
	while( cin >> c )
	{
		if( c == 1 )
		{
			cout	<< "Which sequence?" << endl;
			cin	>> i;
			cout	<< sequences[i] << endl;
		}
		if( c == 2 )
		{
			cout	<< "Which sequences?" << endl;
			cin	>> i >> j;
			dynamic d1( sequences[i], sequences[j] );
			cout	<< d1 << endl;
		}
		cout	<< "Enter command: 1-print, 2-align." << endl;
	}

	return 0;
}
