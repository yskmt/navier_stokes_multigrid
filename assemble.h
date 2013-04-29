// assemble matrix and vector
#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <omp.h>
#include "utils.h"

using namespace std;


// 3d full weighting stencil
// [[1 2 1; 2 4 2; 1 2 1] [2 4 2; 4 8 4; 2 4 2]; [1 2 1; 2 4 2; 1 2 1]];
cdouble fw_stencil[3][3][3] =
	{ {{1.0/16.0,2.0/16.0,1.0/16.0}, {2.0/16.0,4.0/16.0,2.0/16.0},
	   {1.0/16.0,2.0/16.0,1.0/16.0}},
		{{2.0/16.0,4.0/16.0,2.0/16.0}, {4.0/16.0,8.0/16.0,4.0/16.0},
		 {2.0/16.0,4.0/16.0,2.0/16.0}},
			{{1.0/16.0,2.0/16.0,1.0/16.0}, {2.0/16.0,4.0/16.0,2.0/16.0},
			 {1.0/16.0,2.0/16.0,1.0/16.0}}};


void fd_matrix( double** M,
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

#endif //ASSEMBLE_H
