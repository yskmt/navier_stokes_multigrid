// v-cycle
#ifndef V_CYCLE_H
#define V_CYCLE_H

#include "utils.h"
#include "jacobi.h"
#include "assemble.h"
#include "IO.h"

// 3d full weighting stencil
// [[1 2 1; 2 4 2; 1 2 1] [2 4 2; 4 8 4; 2 4 2]; [1 2 1; 2 4 2; 1 2 1]];
cdouble fw_stencil[3][3][3] =
	{ {{1.0/64.0,2.0/64.0,1.0/64.0}, {2.0/64.0,4.0/64.0,2.0/64.0},
	   {1.0/64.0,2.0/64.0,1.0/64.0}},
		{{2.0/64.0,4.0/64.0,2.0/64.0}, {4.0/64.0,8.0/64.0,4.0/64.0},
		 {2.0/64.0,4.0/64.0,2.0/64.0}},
			{{1.0/64.0,2.0/64.0,1.0/64.0}, {2.0/64.0,4.0/64.0,2.0/64.0},
			 {1.0/64.0,2.0/64.0,1.0/64.0}}};


// multigrid v-cycle
double* v_cycle( uint n_dof, cuint nx, cuint ny, cuint nz,
				 cdouble hx, cdouble hy, cdouble hz,
				 cdouble hx2i, cdouble hy2i, cdouble hz2i,
				 cdouble tol, cuint max_iteration, cuint pre_smooth_iteration,
				 cdouble lx, cdouble ly, cdouble lz,
				 cuint level, cuint max_level,
				 double* F,
				 double& Er,
				 double* Uss, double* Vss, double* Wss,
				 cdouble bcs[][6]
				 );

// 3D full weight restriction
void restriction( double* R, double* R_new, cuint I, cuint J, cuint K,
				  cuint I_new, cuint J_new, cuint K_new );

// map from fine to coarse solution
void coarse_map( double* R, double* R_new,
				 cuint nei[][3][3],
				 cuint i, cuint j, cuint k, cuint t_new );

// 3D trilinear interpolation
void interpolation( double* U, double* U_fine,
					cuint I, cuint J, cuint K,
					cuint I_fine, cuint J_fine, cuint K_fine);

// map from coarse to fine solution
void fine_map( double* U, double* U_new,
			   uint box_old[][2][2],
			   uint box_new[][2][2] );

// 0 level v_cycle
double* v_cycle_0( double* Rp,
				   uint n_dof, cuint nx, cuint ny, cuint nz,
				   cdouble hx, cdouble hy, cdouble hz,
				   cdouble hx2i, cdouble hy2i, cdouble hz2i,
				   cdouble tol, cuint max_iteration, cuint pre_smooth_iteration,
				   cdouble width, cdouble length, cdouble height,
				   cuint level, cuint max_level,
				   double& Er,
				   double* Uss, double* Vss, double* Wss,
				   cdouble bcs[][6]				   
				   );

#endif // V_CYCLE_H

