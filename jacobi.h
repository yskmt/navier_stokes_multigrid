// jacobi method
#ifndef JACOBI_H
#define JACOBI_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <omp.h>

using namespace std;

void jacobi( const double tol, const int max_iteration,
			 const unsigned int n_dof,
			 double* u_new,
			 double* u_old,
			 double** M,
			 double* F,
			 double& E,
			 double* R);

double convergence_check ( double** M,
						   double* U,
						   double* F,
						   double* R,
						   const int n_dof);

#endif //JACOBI_H
