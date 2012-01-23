	/*
File:		dna.h
Title:		Class declaration for class "dna", for use in DNA sequence
		manipulations (e.g. alignments).
Author:		Juan Nunez-Iglesias <jnuneziglesias@hotmail.com>

Description:

1. DATA MEMBERS

	The current version uses an enum type for the basic unit of dna, i.e. a
nucleotide. Simple chars were discarded to provide flexibility to use the
extended nucleotide alphabet (15-letter). The class provides conversion
specifications between chars and nucleotides. It also defines which nucleotides
match which, and with what "strength". The strength represents the probability
that two nucleotides match. For example, the strength of a A-N match is 0.25.

	Class dna has three other static elements (on top of the character-
-nucleotide specifications). One is a simple counter (counts the number of dna
objects initialized). The other two are related to the output of the sequence.
They are the printing mode (currently available modes are: FASTA, Nice, and
Raw), and the line wrap (i.e. the length of a line of sequence in printing).
note that the line length is the number of nucleotides, not characters; thus
the number of characters per line will be higher in Nice mode (but not in FASTA
mode).

	In addition to the static data members used described above, each dna
object has three elements:
		- the name of the sequence (string)
		- the nucleotide sequence (vector<nucleotide>)
		- the length of the sequence (int)

2. FUNCTIONS

	2.1. FRIENDS

	Three friends of class dna are declared in this file (and defined in
the source file).
		- operator>>: reads a sequence from an input stream. The
	sequence must be in FASTA format. Invalid characters are discarded.
	Lowercase characters are converted to uppercase before being validated.
		- operator<<: outputs a sequence to a stream in the current
	print mode and using the current line wrap. If the line wrap is set to
	0 or less, only the sequence name is output.
		- compare(): this function compares two nucleotides, and
	returns their matching strength as a type double.

	2.2. CONSTRUCTORS, DESTRUCTOR, AND ASSIGNMENT

	A default constructor, a copy constructor, a destructor, and an assign-
ment operator are defined for the class. Their use is intuitive. One important
feature of the default constructor is that it checks whether a dna object has
been previously initialized using the initialization flag. If no DNA objects
have been initialized previously in the program execution, the information of
nucleotide matching strengths must be hardcoded at this point. This is done by
three private functions:
		- valid_init(): sets the valid nucleotide characters to their
	corresponding nucleotides. A non-valid character corresponds to X.
		- match_init(): sets the nucleotide equivalences (e.g. S is
	equivalent to G or C) and their matching strengths.
		- print_init(): sets the characters corresponding to each
	nucleotide.

	2.3. "GET" FUNCTIONS

	A number of functions are provided to find out the value of certain
members. These functions are:
		- name(): returns a string containing the name of the sequence.
		- letter(): returns the character corresponding to the desired
	nucleotide of the sequence.
		- length(): returns the length of the sequence.
		- pmode(): returns the current printing mode.
		- wrap(): returns the current wrap.
		- n(): returns the number of objects currently in scope.
		- operator[]: returns the nucleotide corresponding to the
	specified position. Performs range checking.

	2.4. "SET" FUNCTIONS

	Functions provided to set elements of the object:
		- name(): Takes a string as input, and sets the sequence name
	to that string.
		- sequence(): Sets the nucleotide sequence from a string
	argument. Invalid characters are discarded. Lowercase is automatically
	converted to uppercase.
		- pmode(): Sets the printing mode to the specified value.
		- wrap(): Sets the line wrap to the specified value.

3. NOTES

	3.1. NUCLEOTIDES

	The current nucleotide implementation is temporary and will be replaced
by a 4-bit variable. The bitwise approach allows dynamic calculation of the
nucleotide matching strengths, rather than the inelegant and error-prone hard-
coding currently in use.

	*/


#ifndef DNA_H
#define DNA_H

#include <iostream>
#include <vector>
#include <string>

enum printmode { FASTA, NICE, RAW };

enum nucleotide { A, C, G, T, R, Y, K, M, S, W, B, D, H, V, N, X };

const int DNAGRP = 10; // Length of group of characters when printing in "NICE"
const int DNAALPHA = 15;	// Length of the nucleotide alphabet
const int DNANAME = 100;	// the starting length for a name string
const int DNASEQ = 100;		// the starting length for a sequence vector
const int DNAWRAP = 60;		// the default line length for printing
const int DNAASCII = 128;	// the number of ASCII characters

class dna
{
		friend istream& operator>>( istream&, dna& );
		friend ostream& operator<<( ostream&, const dna& );
		friend double compare( const nucleotide&, const nucleotide& );
	public:
		// constructors, destructor, and assignment operator
		dna( void );
		dna( const dna& );
		~dna();
		dna& operator=( const dna& );

		// "get" functions:
		string name( void ) const { return name_; }
		char letter( int i ) const { return nucprint_[ sequence_[i] ]; }
		int length( void ) const { return length_; }
		static printmode pmode( void ) { return pmode_; }
		static int wrap( void ) { return wrap_; }
		static int n( void ) { return n_; }
		nucleotide operator[]( int ) const;

		// "set" functions:
		void name( string n ) { name_ = n; }
		void sequence( string s );
		static void pmode( printmode pm ) { pmode_ = pm; }
		static void wrap( int w ) { wrap_ = w; }

	private:
		// copy and destroy functions for use by copy constructor,
		// destructor, and assignment operator
		void copy( const dna& );
		void destroy( void );

		void valid_init( void );
		void match_init( void );
		void print_init( void );

	private:
		string name_;			// the sequence name
		vector<nucleotide> sequence_;	// the DNA sequence
		int length_;			// the sequence length

		static vector<nucleotide> valid_;
					// valid nucleotide characters
		static vector< vector<double> > matching_;
					// which nucleotides match which, and
					// with what value
		static vector<char> nucprint_;
					// what chars correspond to each nucl.
		static printmode pmode_;// the printing mode (fasta, nice, ...)
		static int wrap_;	// the line wrap length for printing
		static int n_;		// the number of dna objects initialized

		static bool init_;	// whether a dna object has been
					// initialized during program execution
};

#endif
