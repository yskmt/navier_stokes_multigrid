// assemble matrix and vector
#ifndef ADVECTION_H
#define ADVECTION_H

#include "utils.h"

using namespace std;

cuint X_DIR=0;
cuint Y_DIR=1;
cuint Z_DIR=2;
cuint XY_DIR=3;
cuint XZ_DIR=4;
cuint YZ_DIR=5;


// treat nonlinear advection terms
void advection( cuint nx, cuint ny, cuint nz );

// generate grid matrix
// dir: direction of velocity: 0:x 1:y 2:z
void grid_matrix( boost::multi_array<double, 3>& Ue,
				  cuint nx, cuint ny, cuint nz,
				  cuint dir,
				  cdouble u0, // x0
				  cdouble u1, // xl
				  cdouble u2, // y0
				  cdouble u3, // yl
				  cdouble u4, // z0
				  cdouble u5 ); //zl

// average values in whatever direction
void average( const boost::multi_array<double, 3>& Ue,
			  boost::multi_array<double, 3>& Ua,
			  cuint dir );

#endif //ASSEMBLE_H
