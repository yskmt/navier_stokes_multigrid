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
cuint fw_stencil[3][3][3] =
	{ {{1,2,1}, {2,4,2}, {1,2,1}},
		{{2,4,2}, {4,8,4}, {2,4,2}},
			{{1,2,1}, {2,4,2}, {1,2,1}}};


void fd_matrix( double** M,
				cuint I, cuint J, cuint K,
				const double dx2i,
				const double dy2i,
				const double dz2i
				);

void load_vector( double* F,
				  cuint n_dof,
				  cuint I,
				  cuint J,
				  cuint K);

int boundary_conditins( cuint n_dof,
						cuint I,
						cuint J,
						cuint K,
						double** M
						);

// 3D full weight restriction
void restriction( double* R, double* R_new, cuint I, cuint J, cuint K);


void coarse_map( double* R, double* R_new,
				 unsigned int nei[][3][3],
				 cuint i, cuint j, cuint k, cuint t_new );

#endif //ASSEMBLE_H
