// assemble matrix and vector
#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#include "utils.h"

using namespace std;

void fd_matrix( double** M,
				cuint I, cuint J, cuint K,
				const double dx2i,
				const double dy2i,
				const double dz2i,
				cuint n_dof
				);

// 2nd order stencil
void fd_matrix_sparse( 	vector<tuple <uint, uint, double> >& M_sp,
						vector<double>& val,
						vector<uint>& col_ind,
						vector<uint>& row_ptr,
						cuint I, cuint J, cuint K,
						const double dx2i,
						const double dy2i,
						const double dz2i,
						cuint n_dof
						);

void load_vector( double* F,
				  cuint n_dof,
				  cuint I,
				  cuint J,
				  cuint K
				  );

int boundary_conditins( cuint n_dof,
						cuint I,
						cuint J,
						cuint K,
						double** M,
						double* F
						);

// add index and value into a sparse matrix
void sparse_add( vector<tuple<uint, uint, double > >& M,
					cuint i, cuint j, cdouble v);

// merge two sorted arrays
void merge(vector<tuple <uint, uint, double> >& left,
		   vector<tuple <uint, uint, double> >& right,
		   cuint n_left, cuint n_right,
		   vector<tuple <uint, uint, double> >& result,
		   vector<tuple <uint, uint, double> >& tmp
		   );


#endif //ASSEMBLE_H
