// implicit viscosity contributions
#ifndef PRESSURE_H
#define PRESSURE_H

#include "utils.h"
#include "assemble.h"
#include "msort.h"
#include "jacobi.h"
#include "advection.h"

using namespace std;


// compute pressure correction
void pressure( boost::multi_array<double, 3>& U,
			   boost::multi_array<double, 3>& V,
			   boost::multi_array<double, 3>& W,
			   double* P,
			   double* Uss, double* Vss, double* Wss,
			   cuint nx, cuint ny, cuint nz,
			   cdouble bcs[][6],
			   cdouble hx, cdouble hy, cdouble hz,
			   cdouble hx2i, cdouble hy2i, cdouble hz2i,
			   cdouble tol, cuint max_iteration );

// build the load vector of pressure equation
void pressure_rhs( double* F, double* Uss, double* Vss, double* Wss,
				   cuint nx, cuint ny, cuint nz,
				   cdouble bcs[][6],
				   cdouble hx, cdouble hy, cdouble hz );


// build a pressure matrix
void pressure_matrix( vector<tuple <uint, uint, double> >& Lp_sp,
					  vector<double>& val,
					  vector<uint>& col_ind,
					  vector<uint>& row_ptr,
					  cuint nx, cuint ny, cuint nz,
					  const double hx2i,
					  const double hy2i,
					  const double hz2i,
					  cuint n_dof
					  );

// compute corrections from pressure values
void compute_corrections( const boost::multi_array<double, 3>& Pr,
						  boost::multi_array<double, 3>& Pr_x,
						  boost::multi_array<double, 3>& Pr_y,
						  boost::multi_array<double, 3>& Pr_z,
						  cuint nx, cuint ny, cuint nz,
						  cdouble hx, cdouble hy, cdouble hz );

#endif
