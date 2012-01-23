	/*
File:		dna.cpp
Title:		Class definitions for class "dna" (declared in dna.h)
Author:		Juan Nunez-Iglesias <jnuneziglesias@hotmail.com>
Description:	See class declaration for description of friend and member
		functions. See below for details on implementation.
	*/

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cctype>
#include "dna.h"

/******************************************************************************/

// input a name and sequence in FASTA format

istream& operator>>( istream& input, dna& s )
{
	int i = 0;	// general use index int
	char c;		// temporary input char

	// skip characters until the signal '>'
	while( (input >> c) && c != '>' ) { ; }

	// if we have reached end-of-file, return the stream
	if( input.eof() )
	{
		return input;
	}

	// erase the current contents of the objects
	s.name_.erase();
	s.sequence_.clear();
	s.length_ = 0;

	// reserve a minimum amount of starting memory for name and sequence
	s.name_.resize( DNANAME );
	s.sequence_.reserve( DNASEQ );

	// input the name
	while( input.get(c) )
	{
		// a name could be split into several lines by having a signal
		// '>' at the start of consecutive lines
		if( c == '\n' )
		{
			if( input.peek() == '>' )
			{
				input.get(c);
				// check that the size has not been exceeded
				if( i == static_cast<int> (s.name_.size()) )
				{
					s.name_.resize( i*2 );
				}
				s.name_[i] = ' ';
				i++;
				continue;
			}
			// if the character immediately following the newline
			// was not a signal '>', end name input
			else
			{
				break;
			}
		}
		if( i == static_cast<int>( s.name_.size() ) )
		{
			s.name_.resize( i*2 );
		}
		s.name_[i] = c;
		i++;
	}

	// Then input the sequence
	while( input.get(c) )
	{
		// Sequence read-in continues until either end-of-file is
		// reached or a new sequence begins (indicated by signal '>')
		if( c == '\n' )
		{
			if( input.peek() == '>' )
			{
				break;
			}
			continue;
		}

		// dna uses only uppercase characters
		c = toupper( c );

		// check that the character being input corresponds to a valid
		// nucleotide
		if( s.valid_[c] != X )
		{
			s.sequence_.push_back( s.valid_[c] );
		}
	}
	s.length_ = s.sequence_.size();

	// This allows the last sequence to be read in a while-loop's condition
	// statement (e.g. while( cin >> dna1 )). Without the clear() statement
	// the input will return having reached end-of-file even though the
	// read-in was successful. Thus eof checking is performed before name
	// and sequence read-in.
	input.clear();

	return input;
}

/******************************************************************************/

// outputs the sequence in the specified format

ostream& operator<<( ostream& output, const dna& s )
{
	// First output the name with a signal '>' in front
	output	<< ">" << s.name_;

	// if the wrap is set to 0 or less, return
	if( s.wrap_ < 1 )
	{
		output	<< endl;
		return output;
	}

	// printing in NICE
	if( s.pmode_ == NICE )
	{
		for( int i = 0; i < s.length_; i++ )
		{
			// jump to a newline and display the current nucleotide
			// when the line wrap value is reached
			if( i % s.wrap_ == 0 )
			{
				output	<< endl;
				output.width(6);
				output	<< i+1 << " ";
			}
			// space groups of nucleotides (default 10)
			else if( i % DNAGRP == 0 )
			{
				output	<< " ";
			}
			// print the char corresponding to the current
			// nucleotide
			output	<< ( s.nucprint_[ s.sequence_[i] ] );
		}
		output	<< endl;
	}

	// printing in FASTA
	else if( s.pmode_ == FASTA )
	{
		for( int i = 0; i < s.length_; i++ )
		{
			// if the line wrap is reached, endline
			if( i % s.wrap_ == 0 )
			{
				output	<< endl;
			}
			// output the char corresponding to the current
			// nucleotide
			output	<< ( s.nucprint_[ s.sequence_[i] ] );
		}
	}

	// printing in RAW simply outputs the sequence
	else if( s.pmode_ == RAW )
	{
		output	<< endl;
		for( int i = 0; i < s.length_; i++ )
		{
			output	<< ( s.nucprint_[ s.sequence_[i] ] );
		}
	}

	return output;
}
/******************************************************************************/

// compare two nucleotides

double compare( const nucleotide& n1, const nucleotide& n2 )
{
	// verify that a dna object (and consequently nucleotide match
	// strengths) has been initialized
	if ( !dna::init_ )
	{
		dna d;
	}

	// return the match strength of the two nucleotides
	return dna::matching_[n1][n2];
}

/******************************************************************************/

// definition of static elements of class dna

	// vector mapping each valid char to a nucleotide
	vector<nucleotide> dna::valid_( DNAASCII, X );

	// 2D vector defining the match strengths between nucleotides
	vector< vector<double> > dna::matching_( DNAALPHA, DNAALPHA );

	// vector mapping nucleotides to corresponding chars
	vector<char> dna::nucprint_( DNAALPHA );

	// a data type to specify printing mode for DNA sequences
	printmode dna::pmode_( NICE );

	// the number of characters to print in each line when printing
	int dna::wrap_( DNAWRAP );

	// definition of dna sequence counter
	int dna::n_( 0 );

	// definition of the initialization flag
	bool dna::init_( false );

/******************************************************************************/

// default constructor for class dna

dna::dna( void )
	:
	sequence_( 0 ),
	length_( 0 )
{
	// hardcode nucleotide alphabet equivalences if nor initialized before
	if( !init_ )
	{
		valid_init(); // defines which characters are valid and maps
				// them to nucleotides
		match_init(); // defines which nucleotides match which, and
				// with what strength
		print_init(); // defines what characters a nucleotide corres-
				// ponds to
		init_ = true;
	}

	// Increment the counter
	n_++;
}

/******************************************************************************/

// copy constructor for class dna

dna::dna( const dna& d1 )
{
	copy( d1 );
	n_++;
}

/******************************************************************************/

// destructor for class dna

dna::~dna()
{
	destroy();
	n_--;
}

/******************************************************************************/

// assignment operator for class dna

dna& dna::operator=( const dna& d1 )
{
	if( this != &d1 )
	{
		destroy();
		copy( d1 );
	}
	return *this;
}

/******************************************************************************/

// return the specified nucleotide

nucleotide dna::operator[]( int i ) const
{
	// check the nucleotide range
	if( i > length_-1 || i < 0 )
	{
		cerr	<< "ERROR: tried to access nucleotide outside of "
			<< "sequence range.\n Sequence length: "
			<< length_
			<< "; Tried to access:" << i << endl;
		exit( EXIT_FAILURE );
	}
	return sequence_[i];
}

/******************************************************************************/

// set the sequence to a specified string

void dna::sequence( string s )
{
	// clear the current contents of the sequence.
	sequence_.clear();

	// set the length
	length_ = static_cast<int> (s.length());

	// set the sequence
	for( int i = 0; i < length_; i++ )
	{
		// convert lowercase to uppercase
		s[i] = toupper( s[i] );
		// discard invalid characters
		if( valid_[s[i]] != X )
		{
			sequence_.push_back( valid_[s[i]] );
		}
	}

	return;
}

/******************************************************************************/

// copy an input dna sequence to the calling sequence

void dna::copy( const dna& d1 )
{
	name_ = d1.name_;
	sequence_ = d1.sequence_;
	length_ = d1.length_;
}

/******************************************************************************/

// destroy the current sequence. Since all dynamic memory allocation is
// currently handled by other objects (string and vector), nothing needs to be
// done by destroy()

void dna::destroy( void )
{
}

/******************************************************************************/

// hardcode the valid characters for nucleotide representation

void dna::valid_init( void )
{
	valid_['A'] = A;
	valid_['C'] = C;
	valid_['G'] = G;
	valid_['T'] = T;
	valid_['R'] = R;
	valid_['Y'] = Y;
	valid_['K'] = K;
	valid_['M'] = M;
	valid_['S'] = S;
	valid_['W'] = W;
	valid_['B'] = B;
	valid_['D'] = D;
	valid_['H'] = H;
	valid_['V'] = V;
	valid_['N'] = N;
}

/******************************************************************************/

// hardcode the match strengths of different nucleotides. Sorted by decreasing
// match strength, and by nucleotide order (i.e. the number of actual
// nucleotides represented by the symbol)

void dna::match_init( void )
{
	matching_[A][A] = matching_[C][C] = matching_[G][G] = matching_[T][T]
	= 1.0;

	matching_[A][R] = matching_[A][W] = matching_[A][M]
	= matching_[R][A] = matching_[W][A] = matching_[M][A]
	= matching_[C][Y] = matching_[C][S] = matching_[C][M]
	= matching_[Y][C] = matching_[S][C] = matching_[M][C]
	= matching_[G][R] = matching_[G][S] = matching_[G][K]
	= matching_[R][G] = matching_[S][G] = matching_[K][G]
	= matching_[T][Y] = matching_[T][W] = matching_[T][K]
	= matching_[Y][T] = matching_[W][T] = matching_[K][T]
	= 0.5;

	matching_[R][R] = matching_[Y][Y] = matching_[K][K] = matching_[M][M]
	= matching_[S][S] = matching_[W][W]
	= 0.5;

	matching_[A][D] = matching_[A][H] = matching_[A][V]
	= matching_[D][A] = matching_[H][A] = matching_[V][A]
	= matching_[C][B] = matching_[C][H] = matching_[C][V]
	= matching_[B][C] = matching_[H][C] = matching_[V][C]
	= matching_[G][B] = matching_[G][D] = matching_[G][V]
	= matching_[B][G] = matching_[D][G] = matching_[V][G]
	= matching_[T][B] = matching_[T][D] = matching_[T][H]
	= matching_[B][T] = matching_[D][T] = matching_[H][T]
	= 1/3;

	matching_[R][D] = matching_[R][V]
	= matching_[D][R] = matching_[V][R]
	= matching_[Y][B] = matching_[Y][H]
	= matching_[B][Y] = matching_[H][Y]
	= matching_[K][B] = matching_[K][D]
	= matching_[B][K] = matching_[D][K]
	= matching_[M][H] = matching_[M][V]
	= matching_[H][M] = matching_[V][M]
	= matching_[S][B] = matching_[S][V]
	= matching_[B][S] = matching_[V][S]
	= matching_[W][D] = matching_[W][H]
	= matching_[D][W] = matching_[H][W]
	= 1/3;

	matching_[B][B] = matching_[D][D] = matching_[H][H] = matching_[V][V]
	= 1/3;

	matching_[A][N] = matching_[C][N] = matching_[G][N] = matching_[T][N]
	= matching_[N][A] = matching_[N][C] = matching_[N][G] = matching_[N][T]
	= 0.25;

	matching_[R][N] = matching_[Y][N] = matching_[K][N] = matching_[M][N]
	= matching_[N][R] = matching_[N][Y] = matching_[N][K] = matching_[N][M]
	= matching_[S][N] = matching_[W][N]
	= matching_[N][S] = matching_[N][W]
	= 0.25;

	matching_[B][N] = matching_[D][N] = matching_[H][N] = matching_[V][N]
	= matching_[N][B] = matching_[N][D] = matching_[N][H] = matching_[N][V]
	= 0.25;

	matching_[N][N] = 0.25;

	matching_[R][B] = matching_[R][H]
	= matching_[B][R] = matching_[H][R]
	= matching_[Y][D] = matching_[Y][V]
	= matching_[D][Y] = matching_[V][Y]
	= matching_[K][H] = matching_[K][V]
	= matching_[H][K] = matching_[V][K]
	= matching_[M][B] = matching_[M][D]
	= matching_[B][M] = matching_[D][M]
	= matching_[S][D] = matching_[S][H]
	= matching_[D][S] = matching_[H][S]
	= matching_[W][B] = matching_[W][V]
	= matching_[B][W] = matching_[V][W]
	= 1/6;

}

/******************************************************************************/

// hardcode the characters equivalent to each nucleotide

void dna::print_init( void )
{

	nucprint_[A] = 'A';
	nucprint_[C] = 'C';
	nucprint_[G] = 'G';
	nucprint_[T] = 'T';
	nucprint_[R] = 'R';
	nucprint_[Y] = 'Y';
	nucprint_[K] = 'K';
	nucprint_[M] = 'M';
	nucprint_[S] = 'S';
	nucprint_[W] = 'W';
	nucprint_[B] = 'B';
	nucprint_[D] = 'D';
	nucprint_[H] = 'H';
	nucprint_[V] = 'V';
	nucprint_[N] = 'N';

}

/******************************************************************************/

