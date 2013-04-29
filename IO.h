// assemble matrix and vector
#ifndef IO_H
#define IO_H

#include "utils.h"

using namespace std;


// write out the sparse matrix
int write_matrix(cuint P,
				 cuint Q,
				 double** U,
				 char* file_name);

// write out the sparse matrix
int write_vector( cuint P,
				  double* F,
				  char* file_name);

// write out the results
int write_results( double* u,
				   cuint n_dof,
				   cuint I,
				   cuint J,
				   cuint K,
				   cdouble dx,
				   cdouble dy,
				   cdouble dz,
				   cuint level
				   );


#endif //IO_H
