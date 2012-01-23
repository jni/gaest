	/*
File:		gaest.cpp
Title:		EST clustering program based on Genetic Algorithms
Author:		Juan Nunez-Iglesias <jnuneziglesias@hotmail.com>

Description:	This program makes use of MIT's GAlib (the Genetic Algorithm
		library written by Matthew Wall) to cluster EST sequences
		by similarity.

		Type -h for command-line descriptions.

		See end of this file for implementation details.

		Modules needed: dna.h, dna.cpp, dynamic.h, dynamic.cpp. GAlib
		must be installed.

		NOTE: The hash table implementation is provided in <hash_map>,
		which is not currently part of the C++ STL. It is expected to
		be included in the next version though.

	*/

/******************************************************************************/

// imported files, symbolic constants

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <hash_map>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <ga/ga.h>
#include "dna.h"
#include "dynamic.h"

// default values for some program parameters
const float LOAD = 0.5;
const char* PARAMFILE = "gaparam.in";
const int MAXSIZE = 1000;

/******************************************************************************/

// declaration of the genome class

template class GA1DArrayGenome<int>;

// declaration of user-defined functions for use with the GA

float objective( GAGenome& );
void initializer( GAGenome& );
int mutator( GAGenome&, float );

// declaration of "helper" functions

int traversecluster( vector<bool>&, vector< list<int> >&, int,
	bool print = false, ostream& output = cout, bool namesonly = false );
void error( const string&, const string& );
void check( int, int );
void printtime( double, ostream& );

// declaration of global vector<> of dna sequences and a 2D vector (the 2nd
// dimension is a hash table to keep space usage down) containing all
// previously determined edges
vector<dna> sequences;
vector< hash_map<int, bool> > scores;

// declaration of global variables for trace statistics
bool trace = false;
string tracefilename = "gaesttrace.out";
ofstream tracefile;
double numaligned = 0;
double expected = 0;
time_t start, end;
double timediff;

/******************************************************************************/

// main program

int main( int argc, char** argv )
{
	// vector of strings to hold command-line arguments
	vector< string > arguments( argc );

	// transfer the arguments to string format
	for( int i = 0; i < argc; i++ )
	{
		arguments[i] = argv[i];
	}

	// variables that can be modified through the command-line
	float hashload = LOAD;
	int maxsize = MAXSIZE;
	string infile, outfile, paramfile, statsfile;
	bool namesonly = false, stats = false;

	// Strings containing error and help messages
	const string errormsg( "Incorrect option syntax. Use -h for help." );
	const string usage( "\nEST Clustering program. Read in sequences from "
		"stdin in FASTA format, clusters\nthem by similarity and prints"
		" them in clusters to stdout.\n\n"
		"Available options:\n"
		"\t-l(oad) float:\tspecify the expected load of hash tables.\n"
			"\t\t\tfloat must be > 0. Low values use more memory\n"
			"\t\t\tbut they are faster, and vice-versa.\n"
		"\t-s(ize) int:\tspecify the maximum size of the hash tables.\n"
		"\t-stats file:\tprint GA statistics to the specified file.\n"
		"\t-i(nput) file:\tspecify a file from which to read in\n"
			"\t\t\tsequences.\n"
		"\t-o(utput) file:\tspecify a file to which to output the\n"
			"\t\t\tresults.\n"
		"\t-p(arams) file:\tspecify the file from which the GA\n"
			"\t\t\tparameters should be input.\n"
		"\t-n(ames):\tonly output sequence names.\n"
		"\t-t(race) file:\tprint trace statistics to a file.\n"
		"\t-h(elp):\tyou probably know this one already... ;-)\n"
		);


	// process arguments (WILL BE MOVED TO OPTION-HANDLING MODULE)
	for( int i = 1; i < argc; i++ )
	{
		string opt( arguments[i] );

		// change the load as specified on the command line
		if( opt == "-l" || opt == "-load" )
		{
			if( i+1 < argc )
			{
				i++;
				hashload = atof( arguments[i].c_str() );
				if( hashload <= 0 )
				{
					error( arguments[0], errormsg );
				}
			}
			else
			{
				error( arguments[0], errormsg );
			}
			continue;
		}

		// print GA statistics to specified file
		if( opt == "-stats" )
		{
			if( i+1 < argc )
			{
				i++;
				stats = true;
				statsfile = arguments[i];
			}
			else
			{
				error( arguments[0], errormsg );
			}
			continue;
		}

		// set the hash table maximum size to the specified size
		if( opt == "-s" || opt == "-size" )
		{
			if( i+1 < argc )
			{
				i++;
				maxsize = atoi( arguments[i].c_str() );
			}
			else
			{
				error( arguments[0], errormsg );
			}
			continue;
		}

		// set the file from which to read input
		if( opt == "-i" || opt == "-input" )
		{
			if( i+1 < argc )
			{
				i++;	
				infile = arguments[i];
			}
			else
			{
				error( arguments[0], errormsg );
			}
			continue;
		}

		// set the file to which to output the results
		if( opt == "-o" || opt == "-output" )
		{
			if( i+1 < argc )
			{
				i++;
				outfile = arguments[i];
			}
			else
			{
				error( arguments[0], errormsg );
			}
			continue;
		}

		// set the file from which to read GA parameters
		if( opt == "-p" || opt == "-params" )
		{
			if( i+1 < argc )
			{
				i++;
				paramfile = arguments[i];
			}
			else
			{
				error( arguments[0], errormsg );
			}
			continue;
		}

		// set the 'namesonly' flag to output only names
		if( opt == "-n" || opt == "-names" )
		{
			namesonly = true;
			continue;
		}

		// print trace statistics. If no file is specified a default
		// file will be created OR replaced.
		if( opt == "-t" || opt == "-trace" )
		{
			trace = true;
			if( i+1 < argc )
			{
				if( arguments[i+1][0] == '-' )
				{
					tracefile.open( tracefilename.c_str(),
						ios::out );
				}
				else
				{
					i++;
					tracefile.open( arguments[i].c_str(),
						ios::noreplace | ios::out );
					if( !tracefile.good() )
					{
						error( arguments[0], "Could "
						"not open trace file." );
					}
				}
			}
			continue;
		}

		// print a help message
		if( opt == "-h" || opt == "-help" )
		{
			cout	<< arguments[0] << ": " << endl
				<< usage << endl;
			return 0;
		}

		else
		{
			error( arguments[0], errormsg );
		}
	}

	// initialize a dna "buffer" variable for input
	dna temp;

	// read in the sequences from file (if specified) or cin (default)
	if( infile.size() == 0 )
	{
		while( cin >> temp )
		{
			sequences.push_back( temp );
		}
	}
	else
	{
		ifstream inputfile( infile.c_str(), ios::nocreate | ios::in );
		if( !inputfile.good() )
		{
			error( arguments[0], "ERROR: Input file could not be "
				"opened. Program terminated." );
		}
		while( inputfile >> temp )
		{
			sequences.push_back( temp );
		}
		inputfile.close();
	}

	// create n to represent the number of sequences, for clutter-free code
	int n = static_cast<int> (sequences.size());

	if( trace )
	{
		tracefile << "Number of sequences:\t\t" << n << endl;
	}

	// resize the first dimension of the scores table
	scores.resize( n );

	// initialize the genome
	GA1DArrayGenome<int> genome( n, objective );
	genome.initializer( ::initializer );
	genome.mutator( ::mutator );

	// initialize the GA
	GASimpleGA ga( genome );

	// set the GA parameters
	if( paramfile.size() == 0 )
	{
		paramfile = PARAMFILE;
	}
	ga.parameters( paramfile.c_str(), gaFalse );

	// calculate the expected number of alignments to be computed.
	// (dependent on genome size, population size, mutation rate, and
	// number of generations). See description for derivation.
	int popSize = ga.populationSize();
	int nGen = ga.nGenerations();
	float pMut = ga.pMutation();
	double done( 0 );

	if( trace )
	{
		tracefile << "Population size:\t\t" << popSize << endl
			<< "Number of generations:\t\t" << nGen << endl
			<< "Mutation rate:\t\t\t" << pMut << "\n" << endl;
	}

	// Round up total (expected) number of gene evaluations (SEE BELOW
	// FOR DERIVATION)
	int tot_gen_eval = static_cast<int>
		(floor( n*popSize + n*popSize*nGen*pMut ) + 1 );

	for( int i = 0; i < tot_gen_eval; i++ )
	{
		done = done + 2 - 2 * done/(n*(n-1));
	}
	expected = done;

	// resize each scores[] hash table to have the specified load with the
	// expected number of alignments, or to use less than the specified
	// space (whichever is smaller)
	int tablesize = static_cast<int> (done/n/hashload);
	if( tablesize > n )
	{
		tablesize = n;
	}
	if( tablesize > maxsize )
	{
		tablesize = maxsize;
	}

	for( int i = 0; i < n; i++ )
	{
		scores[i].resize( tablesize );
	}

	if( trace )
	{
		tracefile << "Expected number of dynamic programming "
			<< "alignments: " << done/2 << endl
			<< "Calculated tablesize: " << tablesize << endl
			<< "Real Tablesize: " << scores[0].bucket_count()
			<< "\n" << endl;
	}

				// run the GA

	if( trace )
	{
		tracefile << "Starting GA...\n" << endl;
		start = time(NULL);
	}

	ga.initialize();

	if( trace )
	{
		tracefile << "Generation:\tTime:\t\tBest Score:\n" << endl;
	}

	for( int i = 0; i < nGen; i++ )
	{
		if( trace )
		{
			end = time(NULL);
			timediff = difftime( end, start );
			tracefile << i << "\t\t";
			printtime( timediff, tracefile );
			tracefile << "\t\t"
				<< ga.statistics().bestIndividual().evaluate()
				<< endl;
		}
		ga.step();
	}

	// print the GA statistics
	if( stats )
	{
		ga.statistics().write( statsfile.c_str() );
	}

	// print out the clustering of the best individual
	const GA1DArrayGenome<int>& best =
		static_cast< const GA1DArrayGenome<int>& >
		( ga.statistics().bestIndividual() );

	// first create a graph of the clusters (adjacency list representation)
	vector<bool> nodes( n, false );
	vector< list<int> > edges( n );

	for( int i = 0; i < n ; i++ )
	{
		int j = best.gene(i);
		if( scores[i][j] )
		{
			edges[i].push_back(j);
			edges[j].push_back(i);
		}
	}

	// then print the sequences in clusters

	// if a file has been specified:
	if( outfile.size() != 0 )
	{
		ofstream output( outfile.c_str(), ios::noreplace | ios::out );

		// make sure the file can be opened
		while( !output.good() )
		{
			cerr	<< "Error: could not open specified file:"
				<< endl
				<< outfile
				<< endl
				<< "for output. Action?"
				<< endl
				<< "  0. Exit.\n"
				<< "  1. Specify new file.\n"
				<< "  2. Overwrite.\n"
				<< endl;
			int response;
			cin	>> response;
			switch( response )
			{
				case 0 :
					error( arguments[0], "Abnormal exit!" );
				case 1 :
					cin	>> outfile;
					output.open( outfile.c_str(),
						ios::noreplace | ios::out );
					break;
				case 2 :
					output.open( outfile.c_str(),
						ios::out );
					break;
				default:
					cerr	<< "Invalid response." << endl;
			}
		}

		// first print sequences in clusters
		for( int i = 0, j = 0; i < n; i++ )
		{
			if( !nodes[i] && edges[i].size() != 0 )
			{
				output	<< "Cluster " << j << endl;
				j++;
				traversecluster( nodes, edges, i, true, output,
					namesonly );
				output	<< endl;
			}
		}

		// then the sequences that have not been clustered
		output	<< "Unclustered sequences:" << endl;
		for( int i = 0; i < n; i++ )
		{
			if( !nodes[i] )
			{
				output << " " << i << ": ";
				if( namesonly )
				{
					output	<< sequences[i].name();
				}
				else
				{
					output	<< sequences[i];
				}
				output	<< endl;
			}
		}
		output << endl;
	}

	// if no file was specified, print to cout
	else
	{
		for( int i = 0, j = 0; i < n; i++ )
		{
			if( !nodes[i] )
			{
				cout	<< "Cluster " << j;
				j++;
				traversecluster( nodes, edges, i, true, cout,
					namesonly );
				cout << endl;
			}
		}
		cout	<< "Unclustered sequences: " << endl;
		for( int i = 0; i < n; i++ )
		{
			if( !nodes[i] )
			{
				cout	<< " " << i << ": ";
				if( namesonly )
				{
					cout	<< sequences[i].name() << endl;
				}
				else
				{
					cout	<< sequences[i] << endl;
				}
			}
		}
	}

	tracefile.close();

	return 0;
}

/******************************************************************************/

float objective( GAGenome& g )
{
	GA1DArrayGenome<int>& genome = static_cast< GA1DArrayGenome<int>& > (g);
	float total_score = 0;
	vector<bool> nodes( genome.length(), false );
	vector< list<int> > edges( genome.length() );

	// first build the clusters
	for( int i = 0; i < genome.length(); i++ )
	{
		int j = genome.gene(i);
		if( scores[i][j] )
		{
			edges[i].push_back(j);
			edges[j].push_back(i);
		}
	}

	// then traverse each to determine their sizes
	for( int i = 0; i < static_cast< int > ( nodes.size() ); i++ )
	{
		if( !nodes[i] )
		{
			int x = traversecluster( nodes, edges, i ) - 1;
			total_score += static_cast< float >( x * x );
		}
	}

	return total_score;
}

/******************************************************************************/

void initializer( GAGenome& g )
{
	GA1DArrayGenome<int>& genome = static_cast< GA1DArrayGenome<int>& > (g);

	// initialize each gene in the GA genome
	for( int i = 0; i < genome.length(); i++ )
	{
		// align to another sequence at random, but not to itself
		int j = GARandomInt( 0, sequences.size()-1 );
		while( i == j )
		{
			j = GARandomInt( 0, sequences.size()-1 );
		}
		genome.gene(i, j);

		// perform alignment only if they haven't been aligned before
		check( i, j );
	}
	return;
}

/******************************************************************************/

int mutator( GAGenome& g, float rate )
{
	GA1DArrayGenome<int>& genome = static_cast< GA1DArrayGenome<int>& > (g);

	// the total number of mutations is the number of genes times pMut
	int total_mutations = static_cast< int >
		( floor( rate * static_cast< float > (genome.length()) ) );

	// if the total mutations is 0, flip a coin to determine whether a
	// mutation will take place
	if( total_mutations == 0 )
	{
		if( GAFlipCoin( rate * static_cast< float >(genome.length()) ) )
		{
			int i = GARandomInt( 0, genome.length()-1 );
			int j = GARandomInt( 0, genome.length()-1 );
			while( i == j )
			{
				j = GARandomInt( 0, genome.length()-1 );
			}
			check( i, j );
			genome.gene(i, j);
			return 1;
		}
		return 0;
	}

	// otherwise randomly mutate genes until the total number of mutations
	// is reached
	for( int c = 0; c < total_mutations; c++ )
	{
		int i = GARandomInt( 0, genome.length()-1 );
		int j = GARandomInt( 0, genome.length()-1 );
		while( i == j )
		{
			j = GARandomInt( 0, genome.length()-1 );
		}
		check( i, j );
		genome.gene(i, j);
	}
	return total_mutations;
}

/******************************************************************************/

// check whether two sequences have been aligned, and if not, align them

void check( int i, int j )
{
	if( scores[i].find(j) == scores[i].end() )
	{
		dynamic d1( sequences[i], sequences[j], true );
		scores[i][j] = scores[j][i] = d1.significant();
	}
	return;
}

/******************************************************************************/

// traverse the sequence clusters using DFS graph traversal. The function
//	returns the number of nodes traversed. If called without the fourth and
//	fifth arguments, no output will be produced.

int traversecluster( vector<bool>& nodes, vector< list<int> >& edges, int i,
	bool print, ostream& output, bool namesonly )
{
	// if the node has already been visited, return
	if( nodes[i] )
	{
		return 0;
	}

	// print if an ostream is available
	if( print )
	{
		output	<< " " << i << ": ";
		if( namesonly )
		{
			output	<< sequences[i].name()
				<< endl;
		}
		else
		{
			output	<< sequences[i] << endl;
		}
	}

	// count the current node and mark as visited
	int count( 1 );
	nodes[i] = true;

	// traverse the nodes accessible from this node
	for( list<int>::iterator pos( edges[i].begin() );
		pos != edges[i].end(); pos++ )
	{
		count += traversecluster( nodes, edges, *pos, print, output,
			namesonly );
	}

	// return the number of traversed nodes
	return count;
}

/******************************************************************************/

// print an error message and exit

void error( const string& progname, const string& errormsg )
{
	cerr	<< progname << ": " << errormsg << endl;
	exit( EXIT_FAILURE );
}

/******************************************************************************/

void printtime( double time, ostream& output )
{
	int hours, minutes, seconds;

	hours = static_cast<int> ( floor( time / 3600 ) );

	if( hours != 0 )
	{
		output	<< hours << "h";
		minutes = static_cast<int> ( floor( time / 60 ) )
			% ( hours * 60 );
		if( minutes != 0 )
		{
			output	<< minutes << "min";
			seconds = static_cast<int> ( floor( time ) )
				% ( hours * 3600 ) % ( minutes * 60 );
			if( seconds != 0 )
			{
				output	<< seconds << "s";
			}
		}
		else
		{
			seconds = static_cast<int> ( floor( time ) )
				% ( hours * 3600 );
			if( seconds != 0 )
			{
				output << seconds << "s";
			}
		}
	}
	else
	{
		minutes = static_cast<int> ( floor( time / 60 ) );

		if( minutes != 0 )
		{
			output	<< minutes << "min";
			seconds = static_cast<int> ( floor ( time ) )
				% ( minutes * 60 );
			if( seconds != 0 )
			{
				output << seconds << "s";
			}
		}
		else
		{
			seconds = static_cast<int> ( floor ( time ) );
			output	<< seconds << "s";
		}
	}

	return;
}

/******************************************************************************/

/*
			INFORMATION ON THE IMPLEMENTATION

A. General overview of the implementation of the GA and the formulation of the
	EST clustering problem as a GA

		The genome is of size n, where n is the number of sequences to
	be clustered. It is an array of ints, and each gene represents an
	alignment between two sequences: the gene index versus the gene value.
	So for example if gene#20 has a value of 6, the gene represents the
	alignment between sequences 20 and 6.

		Whenever the result of an alignment is obtained, it is stored
	in a 'vector' of hash tables ('hash_map') of size n, where each of the
	tables is of size x * exp / n, where
		x is an adjustment factor to control total space usage, and
		exp is the expected number of alignments (see below for
			derivation).

		If another gene is evaluated that needs the same result (e.g.
	gene#6 with value 20), it can access it from the hash table. This is
	preferred over having to perform the alignment again, which is quite
	expensive (O(N*M), where N and M are the lengths of the sequences). By
	contrast, even with a small hash table (which uses linear chaining),
	the chain is unlikely to get as long as N*M, which would be around
	90,000, as ESTs are about 300bp in length. Furthermore stepping through
	a list is faster than calculating dynamic programming scores, and also
	on average we will only need to search halfway through the list.

		The scoring function for the GA was initially to add the scores
	of the alignments together. However this can have problems as two
	highly similar sequences are likely to bring the score up
	disproportionately, and are likely to cluster together by themselves
	rather than form large clusters with other similar sequences.

		Therefore, the GA objective function was changed to look only
	at whether or not two sequences are close enough to be in the same
	cluster. If they were, they contributed 1 to the score, otherwise 0.

		However, this objective also had problems: a large cluster gave
	the same score as two clusters of half the size. Thus it was possible
	to get a high-scoring clustering with only small clusters.

		The final scoring function is therefore to look at cluster
	size. The score increases geometrically with cluster size, so that a
	large cluster scores much higher than two half-sized clusters.
	

B. Derivation of the expected number of dynamic programming alignments to be
	computed.

		Let 'exp' be the this number. 'exp' will depend on two numbers
	(which are not independent):

			1. n, the total number of DNA sequences to be clustered
			2. t, the total number of GA genes to be evaluated

		Keep in mind that a GA gene evaluation does not imply a new
	alignment, as the results of previous alignments are kept in memory.

		n will depend on input. t is the sum of:

			1. the total number of genes evaluated in the first
		generation. In the first generation every gene must be
		evaluated, and the number of genes is the genes per genome
		(i.e. the genome size, equal to n) multiplied by the number of
		genomes (i.e. the population size, referred to hereafter as
		popSize). Thus:

				n * popSize

			2. The total number of genes evaluated in subsequent
		generations. In these generations only mutated genes need to be
		evaluated. We can obtain an estimation of the number of genes
		evaluated here by multiplying the total number of genes
		considered by the mutation rate (pMut, or the probability that
		a gene will be mutated). The total number of genes can be
		obtained by mutiplying the genes per generation (calculated
		above) by the number of generations (nGen). Thus:

				pMut * nGen * n * popSize

		Therefore the total number of evaluated genes t will be:

			t = n * popSize + pMut * nGen * n * popSize
			  = ( popSize + popSize * nGen * pMut ) * n

		To obtain 'exp' we must take into account the probability that
	a given gene being evaluated has not been evaluated before, which we
	will call p(i) for the i'th gene evaluation. p(0) is 1, and then as
	i increases p decreases to 0.

		p(i) will depend on the number of alignments that have been
	performed to that point. Let this number be referred to as done(i).

		done(i) can be defined by a recursive relationship: done(i) is
	done(i-1) added to twice the probability that a new alignment will have 
	been performed, i.e. 2*p(i). (twice because each new alignment actually
	gives the results for two alignments: the performed alignment and also
	the reciprocal alignment) Thus:

			done(0) = 0;
			done(i) = done(i-1) + 2*p(i-1).

		p(i) can then be expressed in terms of done(i) and the total
	number of possible alignments. The probability of a new alignment being
	performed is the number of alignments that have not been done divided
	by the total number of possible alignments. Furthermore, the number of
	alignments that have not been performed is equal to the total number
	of possible alignments minus the number of alignments that *have* been
	performed, i.e. done(i). Finally, the total number of possible
	alignments is n*(n-1) (since a sequence will never be aligned to
	itself). Thus:

			p(i) = { n*(n-1) - done(i) } / (n*(n-1))
			     = 1 - done(i) / (n*(n-1)).

		Therefore:

			done(0) = 0;
			done(i) = done(i-1) + 2 - 2 * done(i-1) / (n*(n-1)).

		And:

			-->  exp = done(t) <--

		So the expected number of alignments can now be calculated
	using a simple recursive rule.

		Note that the number of dynamic alignments called will be half
	this number since only half of the alignments will need to be perfor-
	med.

									      */

