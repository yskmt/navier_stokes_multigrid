// openmp merge sort
#ifndef MSORT_H
#define MSORT_H

#include "utils.h"

using namespace std;

// comparison function for sorting pairs
int comp_tuples( const tuple<uint, uint, double>& i,
				 const tuple<uint, uint, double>& j );

// merge two sorted arrays
void merge(tuple <uint, uint, double>* left,
		   tuple <uint, uint, double>* right,
		   const int n_left, const int n_right,
		   tuple <uint, uint, double>* result,
		   tuple <uint, uint, double>* tmp );

// mergesort with OpenMP parallelism
void mergesort(tuple <uint, uint, double>* vec,
			   const int threads,
			   const int n,
			   tuple <uint, uint, double>* tmp
			   );


#endif // MSORT_H
