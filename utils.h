// utility functions
#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <omp.h>
#include <tuple>

const double pi=3.14159265359;

using namespace std;

typedef const unsigned int cuint;
typedef unsigned int uint;
typedef const double cdouble;
typedef const int cint;

// comparison function for sorting pairs
int comp_pairs( const tuple<uint, uint, double>& i,
				 const tuple<uint, uint, double>& j );

void three_d_to_one_d( const unsigned int i,
					  const unsigned int j,
					  const unsigned int k,
					  const unsigned int I,
					  const unsigned int J,
					   unsigned int& t );

void one_d_to_three_d( const unsigned int t,
					   const unsigned int I,
					   const unsigned int J,
					   unsigned int& i,
					   unsigned int& j,
					   unsigned int& k);

void get_neighbor( unsigned int t[][3][3],
				   cuint i, cuint j, cuint k,
				   cuint I, cuint J, cuint K );

// get node numbers in a box
void get_box( uint t[][2][2],
				   cuint i, cuint j, cuint k,
			  cuint I, cuint J, cuint K );


#endif //UTILS_H
