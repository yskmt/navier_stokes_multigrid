// assemble matrix and vector
#ifndef IO_H
#define IO_H

#include "utils.h"
#include "advection.h"

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

// write out the results
int write_results( boost::multi_array<double, 3>& U,
 boost::multi_array<double, 3>& V,
 boost::multi_array<double, 3>& W,
 double* P,
					cuint n_dof,
					cuint nx,
					cuint ny,
					cuint nz,
					cdouble hx,
					cdouble hy,
					cdouble hz,
					cuint ts,
					cdouble bcs[][6]
					);


// write out matrix for debuggin purpose
int write_3d_data( boost::multi_array<double, 3>& U,
				   char* file_name );

#endif //IO_H
