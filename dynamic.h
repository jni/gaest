	/*
File:		dynamic.h
Title:		Class declaration for class "dynamic", a dynamic programming
		algorithm for sequence alignment.
Author:		Juan Nunez-Iglesias <jnuneziglesias@hotmail.com>

Description:

1. DATA MEMBERS

	The dynamic object contains a large number of data members. The most
central are two dna pointers, which point to the dna sequences being aligned by
the dynamic object. Other data members are classified below.

	1.1. DATA MEMBERS INVOLVED IN THE ALIGNMENT

	The data members used to align the two sequences are:
		- score: the score of the alignment.
		- xlen/ylen: the lengths of the dna sequences being aligned.
	These correspond to the dimensions of the matrices used for the
	dynamic programming algorithm.
		- scrMatrix: the matrix of scores used by the dynamic
	programming algorithm.
		- ptrMatrix: the matrix of pointers used by the dynamic
	programming algorithm.

		- match, msmatch, gapopen, gapxtnd: the rewards and penalties
	used by the algorithm.
		- aligned: flag to indicate whether the sequences have been
	aligned.
		- significance: the number of consecutive matching nucleotides
	needed to achieve significance. The algorithm has the option to stop
	the alignment once significance has been reached.

	1.2. DATA MEMBERS INVOLVED IN PRINTOUT

	Printout of the alignment required traceback through the matrices from
the highest scoring point to the start of the alignment, i.e. where the pointer
is NULL. To this purpose several members are defined:
		- xbegin, ybegin, xend, yend: the coordinates of the start of
	the aligned region and the end of the aligned region, respectively;
		- pathlength: the length of the aligned regions (accounting for
	gaps in the alignment)

		- top, bottom, aligned: the strings containing the sequences
	after being edited for the alignment (i.e. with dashes where gaps are)
	and the string corresponding to the matches ('|' for identity, ':' for
	a partial match).

		- wrap: the line wrap value for printout of the alignment.

2. FUNCTIONS

	2.1. FRIENDS

	The output operator<< is defined for class dynamic. The output is
performed from the output strings top, bottom and aligned. The wrap member
specifes line wrap value. The output prints out the names of the aligned
sequences, the score of the alignment, and then the alignment. A warning is
displayed if the alignment is not significant.

	2.2. CONSTRUCTORS, DESTRUCTOR AND ASSIGNMENT

	Three constructors are defined for class dynamic:
		- default constructor: initializes all data members that are
	independent of the sequences.
		- initialization constructor: requires the two sequences as
	arguments and initializes all the variables needed for alignment, and
	aligns the sequences. The printout variables are not computed since
	only a few applications will need printing of the alignment.
		- copy constructor: All items are copied.

	The destructor and assignment operator for the class work as would be
expected.

	2.3. "GET" FUNCTIONS

	dna1() and dna2() return references to the dna elements aligned in the
object. No checking is performed to verify that the pointers are not NULL.
(therefore, usage of the initialization constructor, and not the default
constructor, is recommended wherever possible).

	score() returns the score of the alignment.

	aligned() verifies that the sequences have been aligned.

	significant() verifies that the alignment is significant (i.e. the
sequences are more similar than a threshold). The threshold value for
significance is defined below, and can be changed by the user as an optional
argument to the constructor.

	2.4. "SET" FUNCTIONS

	input() allows the user to change the sequences being held by the
object. This function is provided for circumstances in which the default
constructor for the class had to be called.

	match(), msmatch(), gapopen() and gapxtnd() modify the alignment
rewards and penalties as specified. Similarly for wrap().

	significance() allows the user to change the significance level of the
alignment. The significance level is defined as the number of consecutive
matching nucleotides needed in the two sequences to achieve significance.

	2.5. OTHER FUNCTIONS

	The function tracepath() traces the alignment from the highest-scoring
cell in the score matrix to the starting cell. The traces are used for the
output of the alignment.

3. NOTES

	3.1. THE align() FUNCTION

	align() is the core of the dynamic class and consists of the dynamic
programming algorithm for local sequence alignment.

	3.2. FUTURE MODIFICATIONS

	Certain specialized uses of dynamic may only need the score of the
alignment. In these cases only two columns of the score and pointer matrices
are needed at a time. This reduces space necessity considerably. Thus a
different dynamic class may be created for these cases. Ideally a base class
will be created containing the base components of a dynamic programming
algorithm, and specialized classes will be derived from that.

	*/

#ifndef DYNAMIC_H
#define DYNAMIC_H

#include <iostream>
#include <vector>
#include <string>
#include "dna.h"

// default rewards and penalties for alignment
const float DYNMATCH = 1.0;
const float DYNMSMATCH = -2.0;
const float DYNGAPOPEN = -6.0;
const float DYNGAPXTND = -0.2;

const int DYNSIG = 40;	// default minimal length of aligned region considered
			// significant

const int DYNWRAP = 60;	// default printing wrap value

// ints to represent pointers in dynamic programming algorithm
const int PTRNULL = 0;
const int PTRLEFT = 1;
const int PTRUP = 2;
const int PTRDIAG = 3;

class dynamic
{
		friend ostream& operator<<( ostream&, dynamic& );
	public:
		// constructors, destructor, assignment operator
		dynamic( float m = DYNMATCH, float mm = DYNMSMATCH,
			float go = DYNGAPOPEN, float gx = DYNGAPXTND,
			int sl = DYNSIG );
		dynamic( dna&, dna&, bool s = false, int sl = DYNSIG,
			float m = DYNMATCH, float mm = DYNMSMATCH,
			float go = DYNGAPOPEN, float gx = DYNGAPXTND );
		dynamic( const dynamic& );
		dynamic& operator=( const dynamic& );
		~dynamic();

		// "get" functions
		dna& dna1( void ) const { return (*dna1ptr_); }
		dna& dna2( void ) const { return (*dna2ptr_); }
		float score( void ) const { return score_; }
		int pathlength() const { return pathlength_; }
		bool aligned( void ) const { return aligned_; }
		bool significant( void );

		// "set" functions
		void input( dna&, dna&, bool s = false );
		void match( float m ) { match_ = m; aligned_ = false; }
		void msmatch( float m ) { msmatch_ = m; aligned_ = false; }
		void gapopen( float g ) { gapopen_ = g; aligned_ = false; }
		void gapxtnd( float g ) { gapxtnd_ = g; aligned_ = false; }
		void wrap( int w ) { wrap_ = w; }
		void significance( int n ) { significance_ = n; }

		// other functions
		void tracepath( void );
	private:
		// alignment function for use on initialization
		void align( bool s = false );

		// copy and destroy functions used by copy constructor,
		// assignment operator and destructor.
		void copy( const dynamic& );
		void destroy( void );

		// helper function for align()
		inline int max( float*, int );

	private:
		dna* dna1ptr_;	// the first (x-value) dna sequence
		dna* dna2ptr_;	// the second (y-value) dna sequence

		float score_;		// the final score of the alignment

		int xlen_, ylen_;	// length of _dna1 and _dna2 sequences

		vector< vector<float> >
			scrMatrix_;	// the score matrix
		vector< vector<int> >
			ptrMatrix_;	// the pointer matrix

		int xbegin_, ybegin_;	// coordinates of start of alignment
		int xend_, yend_;	// coordinates of end of alignment
		int pathlength_;	// the length of the aligned region

		string top_;		// the aligned region of _dna1
		string bottom_;		// the aligned region of _dna2
		string align_;		// match sequence of aligned region
		int wrap_;		// line wrap value to print alignment

		float match_;		// the alignment rewards and penalties
		float msmatch_;		// (positive value = reward)
		float gapopen_;		// (negative value = penalty)
		float gapxtnd_;

		bool aligned_;		// have the dna sequences been aligned?
		int significance_;	// what is the minimum length of a
					// significant alignment?
};

#endif
