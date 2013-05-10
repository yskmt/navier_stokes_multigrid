// advection contributions
#ifndef ADVECTION_H
#define ADVECTION_H

#include "utils.h"

using namespace std;

// set up initial conditions
void initial_conditions( boost::multi_array<double, 3>& U,
						 boost::multi_array<double, 3>& V,
						 boost::multi_array<double, 3>& W,
						 double* P,
						 cuint nx, cuint ny, cuint nz );

// treat nonlinear advection terms
void advection( boost::multi_array<double, 3>& U,
				boost::multi_array<double, 3>& V,
				boost::multi_array<double, 3>& W,
				cuint nx, cuint ny, cuint nz,
				cdouble hx, cdouble hy, cdouble hz,
				cdouble dt,
				cdouble bcs[][6]
				);

// generate grid matrix
// dir: direction of velocity: 0:x 1:y 2:z
void grid_matrix( boost::multi_array<double, 3>* U,
				  boost::multi_array<double, 3>* Ue,
				  cuint nx, cuint ny, cuint nz,
				  cuint dir,
				  cdouble* bc ); // x0 xl y0 yl z0 zl

// average values in whatever direction
void average( const boost::multi_array<double, 3>& Ue,
			  boost::multi_array<double, 3>& Ua,
			  cuint dir );

// get maximum value of the 3d array
double max_3d_array( const boost::multi_array<double, 3>& U );


// get upwinding differences
void upwind_difference( const boost::multi_array<double, 3>& U,
						boost::multi_array<double, 3>& Ud,
						cuint dir );

// get 1d staggered difference from cell vertices into cell center
void staggered_first_difference( const boost::multi_array<double, 3>& UV,
								 boost::multi_array<double, 3>& UV_x,
								 cdouble h,
								 cuint dir );

// get central first difference at center of element
void central_first_difference( const boost::multi_array<double, 3>& U2,
							   boost::multi_array<double, 3>& U2_x,
							   cdouble h,
							   cuint dir );

// get mixed edge values
void calculate_edge_values( const boost::multi_array<double, 3>& Ue,
							const boost::multi_array<double, 3>& Ve,
							const boost::multi_array<double, 3>& We,
							boost::multi_array<double, 3>& UV,
							boost::multi_array<double, 3>& UW,
							boost::multi_array<double, 3>& VW,
							cuint nx, cuint ny, cuint nz);

// consolidate advection terms
void consolidate_advection( boost::multi_array<double, 3>& U,
							boost::multi_array<double, 3>& V,
							boost::multi_array<double, 3>& W,
							boost::multi_array<double, 3>& U2_x,
							boost::multi_array<double, 3>& V2_y,
							boost::multi_array<double, 3>& W2_z,
							boost::multi_array<double, 3>& UV_y,
							boost::multi_array<double, 3>& UW_z,
							boost::multi_array<double, 3>& VU_x,
							boost::multi_array<double, 3>& VW_z,
							boost::multi_array<double, 3>& WU_x,
							boost::multi_array<double, 3>& WV_y,
							cuint nx, cuint ny, cuint nz,
							cdouble dt );

#endif //ASSEMBLE_H
