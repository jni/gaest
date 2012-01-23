	/*
File:		dynamic.cpp
Title:		Definition of class "dynamic", dynamic programming algorithm for
		sequence alignment.
Author:		Juan Nunez-Iglesias <jnuneziglesias@hotmail.com>

Description:	The dynamic programming algorithm is O(n*m) and implemented in
		the function align(). align() has the option of stopping as
		soon as a high enough score to be considered significant is
		reached.
		The default alignment simply marks the end coordinates of the
		alignment, which is not sufficient for printing. Thus if
		printing is necessary the function tracepath() can be called.
		Note that printing is safe, i.e. operator<< will call
		tracepath() if it hasn't been done manually.
		For further details on individual functions see the function
		definition.
	*/


#include "dna.h"
#include "dynamic.h"
#include <vector>
#include <iostream>
#include <cstdlib>

/******************************************************************************/

// simple comparison operators ( <, >, ==, <=, >= ) for class dynamic, compare
// the scores of the two "dynamic" objects

inline bool operator<( const dynamic& d1, const dynamic& d2 )
	{ return ( d1.score() < d2.score() ); }
inline bool operator>( const dynamic& d1, const dynamic& d2 )
	{ return ( d1.score() > d2.score() ); }
inline bool operator==( const dynamic& d1, const dynamic& d2 )
	{ return ( d1.score() == d2.score() ); }
inline bool operator<=( const dynamic& d1, const dynamic& d2 )
	{ return ( d1.score() <= d2.score() ); }
inline bool operator>=( const dynamic& d1, const dynamic& d2 )
	{ return ( d1.score() >= d2.score() ); }

/******************************************************************************/

// outputs the aligned regions of the DNA sequences

ostream& operator<<( ostream& output, dynamic& d1 )
{
	// check that the sequences have been aligned
	if( !d1.aligned() )
	{
		// if the pointers are not NULL, align the sequences
		if( d1.dna1ptr_ != 0 && d1.dna2ptr_ != 0 )
		{
			d1.align();
		}
		// otherwise stop output, printing error message
		else
		{
			cerr	<< "Tried to output uninitialized 'dynamic'..."
				<< endl;
			return output;
		}
	}

	// print a message if the alignment is not significant
	if( !d1.significant() )
	{
		output	<< "WARNING: The alignment is not considered "
			<< "significant." << endl;
	}

	// create dna refs from dna pointers
	dna& dna1( *d1.dna1ptr_ ), & dna2( *d1.dna2ptr_ );

	// output names and score of the alignment
	output	<< "Top sequence: " << dna1.name() << endl
		<< "Bottom sequence: " << dna2.name() << endl
		<< "Score: " << d1.score_ << endl;

	// if the linelength is 0 or less, do not print the aligned regions
	if( d1.wrap_ < 1 )
	{
		return output;
	}

	// trace the alignment if it has not been done already
	if( !d1.pathlength_ )
	{
		d1.tracepath();
	}

	// output the aligned regions 'wrap' characters at a time
	for( int i = 0; i < d1.pathlength_; i += d1.wrap_ )
	{
		output.width(6);
		output	<< d1.xbegin_ + i + 1 << "  "	// top sequence index
			<< d1.top_.substr( i, d1.wrap_ )
			<< endl;
		output.width(6);
		output	<< i + 1 << "  "		// alignment index
			<< d1.align_.substr( i, d1.wrap_ )
			<< endl;
		output.width(6);
		output	<< d1.ybegin_ + i + 1 << "  "	// bottom seq. index
			<< d1.bottom_.substr( i, d1.wrap_ )
			<< endl;
		output	<< endl;
	}
	return output;
}

/******************************************************************************/

// default constructor for class dynamic (AVOID USE IN NORMAL CIRCUMSTANCES)
//	(initializing constructor is safer and more efficient)
//	provided solely for unavoidable circumstances, such as when an array
//	of 'dynamic' objects is needed

dynamic::dynamic( float m, float mm, float go, float gx, int sl )
	:
	// initialize dna pointers to 0
	dna1ptr_( 0 ), dna2ptr_( 0 ),

	// initialize score to 0
	score_( 0 ),

	// initialize matrix dimensions to 0
	xlen_( 0 ), ylen_( 0 ),

	// initialize matrices to 0
	scrMatrix_( 0 ), ptrMatrix_( 0 ),

	// initialize all coordinates to 0
	xbegin_( 0 ), ybegin_( 0 ),
	xend_( 0 ), yend_( 0 ),

	// initialize path length to 0
	pathlength_( 0 ),

	// initialize printing wrap value ( aka linelength )
	wrap_( DYNWRAP ),

	// initialize penalties and rewards for alignment
	match_( m ), msmatch_( mm ), gapopen_( go ), gapxtnd_( gx ),

	// set the aligned flag to false
	aligned_( false ),

	// set the significance level
	significance_( sl )
{
}
	
/******************************************************************************/

// initializing constructor for class dynamic

dynamic::dynamic( dna& d1, dna& d2, bool s, int sl,
	float m, float mm, float go, float gx )
	:
	// initialize dna elements
	dna1ptr_( &d1 ), dna2ptr_( &d2 ),

	// initialize score to 0
	score_( 0 ),

	// initialize x and y lengths for matrices
	xlen_( dna1ptr_->length() ), ylen_( dna2ptr_->length() ),

	// initialize matrices to appropriate sizes
	scrMatrix_( xlen_, ylen_ ), ptrMatrix_( xlen_, ylen_ ),

	// initialize all other coordinates to 0
	xbegin_( 0 ), ybegin_( 0 ), xend_( 0 ), yend_( 0 ),

	// initialize path length to 0
	pathlength_( 0 ),

	// initialize printing wrap value ( aka linelength )
	wrap_( DYNWRAP ),

	// initialize penalties and rewards for alignment
	match_( m ), msmatch_( mm ), gapopen_( go ), gapxtnd_( gx ),

	// set the significance level
	significance_( sl )
{
	// align the sequences
	align( s );
}

/******************************************************************************/

// copy constructor for class dynamic

dynamic::dynamic( const dynamic& d1 )
{
	copy( d1 );
}

/******************************************************************************/

// assignment operator for class dynamic

dynamic& dynamic::operator=( const dynamic& d1 )
{
	if( *this == d1 ) // avoid self-assignment
	{
		return *this;
	}
	destroy();
	copy( d1 );
	return *this;
}
/******************************************************************************/

// destructor for class dynamic

dynamic::~dynamic()
{
	destroy();
}

/******************************************************************************/

// input two dna sequences to align.

void dynamic::input( dna& d1, dna& d2, bool s )
{
	// set the dna pointers
	dna1ptr_ = &d1;
	dna2ptr_ = &d2;

	// set the x and y lengths
	xlen_ = d1.length();
	ylen_ = d2.length();

	// make sure any previous matrices are properly deleted
	destroy();

	// initialize score and pointer matrices

	// initialize first dimension
	scrMatrix_.resize( xlen_ );
	ptrMatrix_.resize( xlen_ );

	// initialize second dimension
	for( int i = 0; i < xlen_; i++ )
	{
		scrMatrix_[i].resize( ylen_ );
		ptrMatrix_[i].resize( ylen_ );
	}

	// align the sequences
	align( s );
}

/******************************************************************************/

// align the two dna sequences using the dynamic programming algorithm
//	NOTE: The pointer matrix is implemented as a matrix of ints. This allows
//	selection of array values using the pointers directly, and is also
//	easier to debug. The symbolic constants for NULL, UP, DIAG, and LEFT are
//	defined in the header file.

void dynamic::align( bool s )
{
	// dna refs to avoid pointer notation
	dna& dna1( *dna1ptr_ ), & dna2( *dna2ptr_ );

	// declare index ints for general use
	register int i, j;

	double res = 0.0; // cache variable stores nucleotide comparison
			     // results

	float allscores[4];	// an array containing the possible scores for
				// any given cell in the score matrix (only the
				// maximum will be chosen)

	allscores[PTRNULL] = 0;	// since this is a local alignment, 0 is always
				// a choice for a cell

	// initialize first row of scores and pointers
	for( i = 0; i < xlen_; i++ )
	{
		// the score of the cell will be a match if the nucleotides
		// match, 0 otherwise
		scrMatrix_[i][0] = ( res = compare( dna1[i], dna2[0] ) )?
					res * match_ : 0;

		// the pointers will always be NULL
		ptrMatrix_[i][0] = PTRNULL;
	}

	// initialize first column of scores and pointers
	for( j = 1; j < ylen_; j++ )
	{
		// [ as above ]
		scrMatrix_[0][j] = ( res = compare( dna1[0], dna2[j] ) )?
					res * match_ : 0;
		ptrMatrix_[0][j] = PTRNULL;
	}

	// calculate the local alignment score of each cell in the matrix,
	// updating the pointers along the way.
	for( j = 1; j < ylen_; j++ )
	{
		for( i = 1; i < xlen_; i++ )
		{
			// the score coming from the LEFT will be the score
			// from the LEFT adjacent cell, plus a gap opening
			// penalty if a new gap is being formed, or a gap
			// extension penalty otherwise
			allscores[PTRLEFT] = scrMatrix_[i-1][j] +
					( ( ptrMatrix_[i-1][j] == PTRLEFT ) ?
					gapxtnd_ : gapopen_ );

			// [ as above, replace LEFT with UP ]
			allscores[PTRUP] = scrMatrix_[i][j-1] +
					( ( ptrMatrix_[i][j-1] == PTRUP ) ?
					gapxtnd_ : gapopen_ );

			// the score coming from the diagonal will be the score
			// from the diagonally adjacent cell plus a match if
			// the nucleotides match, a mismatch otherwise
			allscores[PTRDIAG] = scrMatrix_[i-1][j-1] + (
					( res = compare( dna1[i], dna2[j] ) ) ?
					( res * match_ ) : msmatch_ );

			// the pointer indicates the origin of the highest of
			// the four scores
			ptrMatrix_[i][j] = max( allscores, 4 );

			// the score will be the highest of the four scores, as
			// indicated by the pointer
			scrMatrix_[i][j] = allscores[ ptrMatrix_[i][j] ];

			// substitute the maximum score if it has been surpassed
			if( scrMatrix_[i][j] > score_ )
			{
				score_ = scrMatrix_[i][j];
				xend_ = i;
				yend_ = j;

				// if we are only interested in whether the
				// alignment is significant, (that's what s
				// is for... See header documentation), and
				// we have reached significance, exit.
				if( s && significant() )
				{
					return;
				}
			}
		}
	}

	// set the aligned flag to true
	aligned_ = true;

	return;
}

/******************************************************************************/

// copy the contents of a dynamic object to this one (used by copy constructor
//	and assignment operator).

void dynamic::copy( const dynamic& d1 )
{
	// copy elements that must be initialized
	pathlength_ = d1.pathlength_;
	wrap_ = d1.wrap_;
	match_ = d1.match_;
	msmatch_ = d1.msmatch_;
	gapopen_ = d1.gapopen_;
	gapxtnd_ = d1.gapxtnd_;

	// copy elements that will be present if sequences have been aligned
	if( d1.aligned_ )
	{
		aligned_ = d1.aligned_;
		significance_ = d1.significance_;

		dna1ptr_ = d1.dna1ptr_;
		dna2ptr_ = d1.dna2ptr_;
		score_ = d1.score_;
		xlen_ = d1.xlen_;
		ylen_ = d1.ylen_;

		// create space for dynamically allocated elements

		scrMatrix_.resize( xlen_ );
		ptrMatrix_.resize( xlen_ );

		for( int i = 0; i < xlen_; i++ )
		{
			scrMatrix_[i].resize( ylen_ );
			ptrMatrix_[i].resize( ylen_ );
	
			// copy dynamically alocated elements

			for( int j = 0; j < ylen_; j++ )
			{
				scrMatrix_[i][j] = d1.scrMatrix_[i][j];
				ptrMatrix_[i][j] = d1.ptrMatrix_[i][j];
			}
		}
		// check whether the path has been traced, and if so copy
		// 	appropriate elements
		if( pathlength_ )
		{
			xbegin_ = d1.xbegin_;
			ybegin_ = d1.ybegin_;
			xend_ = d1.xend_;
			yend_ = d1.yend_;

			top_ = d1.top_;
			bottom_ = d1.bottom_;
			align_ = d1.align_;
		}
	}
}

/******************************************************************************/

// deletes dynamically allocated elements (used by assignment operator and by
//	class destructor).

void dynamic::destroy( void )
{
	scrMatrix_.clear();
	ptrMatrix_.clear();
}

/******************************************************************************/

// is the alignment significant (i.e. are the sequences related?). The function
//	allows for a 5% mismatch if the minimum length region is aligned

bool dynamic::significant( void )
{
	if( !aligned_ )
	{
		return false;
	}
	if( score_ < (significance_ * (match_ + 0.05 * msmatch_)) )
	{
		return false;
	}
	return true;
}

/******************************************************************************/

// tracepath() traces the aligned portions of the sequence. It then stores the
//	regions of the two sequences in strings, as well as an "align" string
//	that has a '|' at every match on the alignment. The complexity is O(L),
//	where L is the length of the aligned region.

void dynamic::tracepath( void )
{
	dna& dna1( *dna1ptr_ ), & dna2( *dna2ptr_ );
	int i = xend_, j = yend_;

	// perform first pass to determine the starting point of the alignment
	while( ptrMatrix_[i][j] != PTRNULL )
	{
		switch( ptrMatrix_[i][j] )
		{
			case PTRDIAG:
				i--;
				j--;
				pathlength_++;
				break;
			case PTRLEFT:
				i--;
				pathlength_++;
				break;
			case PTRUP:
				j--;
				pathlength_++;
				break;
			default:
				cerr	<< "ERROR: Unrecognized pointer value "
					<< "in dynamic programming ptrMatrix.\n"
					<< "x = " << i << "; y = " << j << endl;
				exit( EXIT_FAILURE );
		}
	}

	// set the beginning coordinates
	xbegin_ = i;
	ybegin_ = j;

	// resize the alignment strings to the pathlength
	top_.resize( pathlength_ );
	bottom_.resize( pathlength_ );
	align_.resize( pathlength_ );

	i = xend_;
	j = yend_;

	double res = 0; // cache variable to store computed results

	// then perform second pass to copy sequences into strings
	for( int k = pathlength_-1; k >= 0 ; k-- )
	{
		switch( ptrMatrix_[i][j] )
		{
			case PTRDIAG:
				top_.at(k) = dna1.letter(i);
				bottom_.at(k) = dna2.letter(j);
				res = compare( dna1[i], dna2[j] );
				if( res == 1 )
				{
					align_.at(k) = '|';
				}
				else if ( res == 0 )
				{
					align_.at(k) = ' ';
				}
				else
				{
					align_.at(k) = ':';
				}
				i--;
				j--;
				break;
			case PTRLEFT:
				top_.at(k) = dna1.letter(i);
				bottom_.at(k) = '-';
				align_.at(k) = ' ';
				i--;
				break;
			case PTRUP:
				top_.at(k) = '-';
				bottom_.at(k) = dna2.letter(j);
				align_.at(k) = ' ';
				j--;
				break;

			// the last nucleotide has been reached
			default:
				top_.at(k) = dna1.letter(i);
				bottom_.at(k) = dna2.letter(j);
				res = compare( dna1[i], dna2[j] );
				if( res == 1 )
				{
					align_.at(k) = '|';
				}
				else if( res == 0 )
				{
					align_.at(k) = ' ';
				}
				else
				{
					align_.at(k) = ':';
				}
		}
	}
	return;
}

/******************************************************************************/

// retrieve the maximum value from a float array of length l (returns its index)

int dynamic::max( float* farray, int l )
{
	int imax = 0;
	for( int i = 1; i < l; i++ )
	{
		if( farray[i] > farray[imax] )
		{
			imax = i;
		}
	}
	return imax;
}

/******************************************************************************/
