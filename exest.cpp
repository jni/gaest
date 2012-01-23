

#include <iostream>
#include <vector>
#include <list>
#include <ctime>
#include "dna.h"
#include "dynamic.h"

const int NUMSEQS = 20;

int traversecluster( vector<bool>&, vector< list<int> >&, int );

vector<dna> sequences;
vector< vector<bool> > edges;

int main( void )
{
	dna temp;
	time_t start, end;

	sequences.reserve( NUMSEQS );

	while( cin >> temp )
	{
		sequences.push_back( temp );
	}

	int n = static_cast<int> (sequences.size());

	cout	<< "Number of sequences: " << n << "\n\n" << endl;

	edges.resize( n );
	for( int i = 0; i < n; i++ )
	{
		edges[i].resize( n );
	}

	start = time( NULL );

	for( int i = 0; i < n; i++ )
	{
		for( int j = i+1; j < n; j++ )
		{
			dynamic d1( sequences[i], sequences[j], true );
			if( d1.significant() )
			{
				edges[i][j] = edges[j][i] = true;
			}
		}
	}

	end = time( NULL );

	int totalsecs = static_cast<int>( difftime( end, start ) );

	vector<bool> nodes( n );

	vector< list<int> > adjedges( n );

	for( int i = 0; i < n; i++ )
	{
		for( int j = 0; j < n; j++ )
		{
			if( edges[i][j] )
			{
				adjedges[i].push_back( j );
			}
		}
	}

	float tempscore( 0 ), totscore( 0 );

	for( int i = 0, j = 0; i < n; i++ )
	{
		if( !nodes[i] && adjedges[i].size() != 0 )
		{
			cout	<< "Cluster " << j << endl << " ";
			tempscore = 0;
			tempscore = traversecluster( nodes, adjedges, i );
			cout	<< endl;
			totscore += ( tempscore-1 ) * ( tempscore-1 );
			j++;
		}
	}

	cout	<< "Singletons: " << endl;
	for( int i = 0; i < n; i++ )
	{
		if( !nodes[i] )
		{
			cout	<< i << " ";
		}
	}
	cout	<< "\n\n SCORE: " << totscore << endl
		<< " TIME: " << totalsecs << endl
		<< " ALIGNMENTS: " << n*(n-1)/2 << endl;

	return 0;
}

int traversecluster( vector<bool>& nodes, vector< list<int> >& adjedges, int i )
{
	if( nodes[i] )
	{
		return 0;
	}

	cout	<< i << " ";

	int count = 1;
	nodes[i] = true;

	for( list<int>::iterator pos( adjedges[i].begin() );
		pos != adjedges[i].end(); pos++ )
	{
		count += traversecluster( nodes, adjedges, *pos );
	}
	return count;
}
