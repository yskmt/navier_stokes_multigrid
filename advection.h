// advection contributions
#ifndef ADVECTION_H
#define ADVECTION_H

#include "utils.h"

using namespace std;

// treat nonlinear advection terms
void advection( double* U,
				double* V,
				double* W,
				cuint nx, cuint ny, cuint nz,
				cdouble hx, cdouble hy, cdouble hz,
				cdouble dt,
				cdouble bcs[][6]
				);
	
// generate grid matrix
// dir: direction of velocity: 0:x 1:y 2:z
void grid_matrix( double* U,
				  double* Ue,
				  cuint nx, cuint ny, cuint nz,
				  cuint nxe, cuint nye, cuint nze,
				  cuint dir,
				  cdouble* bc );

// average values in whatever direction
void average( const double* Ue,
			  double* Ua,
			  cuint nxe, cuint nye, cuint nze,
			  cuint nxa, cuint nya, cuint nza,
			  cuint dir );

// get maximum value of the 3d array
double max_3d_array( const boost::multi_array<double, 3>& U );


// get upwinding differences
void upwind_difference( const boost::multi_array<double, 3>& U,
						boost::multi_array<double, 3>& Ud,
						cuint dir );

// get 1d staggered difference
// Ua is an averaged U value at the cell vertices
// get the difference value at the center of cell
void staggered_first_difference( const double* UV,
								 double* UV_x,
								 cuint nx, cuint ny, cuint nz,
								 cuint nx_x, cuint ny_x, cuint nz_x,
								 cdouble h,
								 cuint dir
								 );

// get central first difference at center of element
void central_first_difference( const boost::multi_array<double, 3>& U2,
							   boost::multi_array<double, 3>& U2_x,
							   cdouble h,
							   cuint dir );


// get mixed edge values
void calculate_edge_values( double* Ue,
							double* Ve,
							double* We,
							double* UV,
							double* UW,
							double* VW,
							cuint nx, cuint ny, cuint nz);

// consolidate advection terms
void consolidate_advection( double* U,
							double* V,
							double* W,
							double* U2_x,
							double* V2_y,
							double* W2_z,
							double* UV_y,
							double* UW_z,
							double* VU_x,
							double* VW_z,
							double* WU_x,
							double* WV_y,
							cuint nx, cuint ny, cuint nz,
							cdouble dt );

#endif //ASSEMBLE_H
