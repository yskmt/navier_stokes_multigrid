// jacobi method
#ifndef JACOBI_H
#define JACOBI_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <omp.h>

using namespace std;

const double pi=3.14159265359;

void jacobi( const double tol, const int max_iteration,
			 const unsigned int n_dof,
			 double* u_new,
			 double* u_old,
			 double** M,
			 double* F,
			 double& E);

#endif //JACOBI_H
