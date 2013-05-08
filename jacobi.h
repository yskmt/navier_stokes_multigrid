// jacobi method
#ifndef JACOBI_H
#define JACOBI_H

#include "utils.h"

using namespace std;

void jacobi( cdouble tol,
			 cuint max_iteration,
			 cuint n_dof,
			 double* u_new,
			 double* u_old,
			 double** M,
			 double* F,
			 double& E,
			 double* R);

// sparse jacobi method
void jacobi_sparse( cdouble tol,
					cuint max_iteration,
					cuint n_dof,
					double* U,
					double* U_tmp,
					const vector<double>& val,
					const vector<uint>& col_ind,
					const vector<uint>& row_ptr,
					double* F,
					double& Er,
					double* R);

double convergence_check ( double** M,
						   double* U,
						   double* F,
						   double* R,
						   cuint n_dof
						   );

double convergence_check_sparse ( const vector<double>& val,
								  const vector<uint>& col_ind,
								  const vector<uint>& row_ptr,
								  double* U,
								  double* F,
								  double* R,
								  cuint n_dof);

#endif //JACOBI_H

