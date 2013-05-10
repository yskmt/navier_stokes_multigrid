// implicit viscosity contributions
#ifndef VISCOSITY_H
#define VISCOSITY_H

#include "utils.h"
#include "assemble.h"
#include "msort.h"
#include "jacobi.h"

using namespace std;


// implicitly solve viscosity
void viscosity(  boost::multi_array<double, 3>& U,
				 boost::multi_array<double, 3>& V,
				 boost::multi_array<double, 3>& W,
				 cuint nx, cuint ny, cuint nz,
				 cdouble hx, cdouble hy, cdouble hz,
				 cdouble hx2i, cdouble hy2i, cdouble hz2i,
				 cdouble dt, cdouble nu,
				 cdouble bcs[][6],
				 cdouble tol, cuint max_iteration );


// implicitly solve viscosity
void viscosity(  boost::multi_array<double, 3>& U,
				 boost::multi_array<double, 3>& V,
				 boost::multi_array<double, 3>& W,
				 double* Uss, double* Vss, double* Wss,
				 cuint nx, cuint ny, cuint nz,
				 cdouble hx, cdouble hy, cdouble hz,
				 cdouble hx2i, cdouble hy2i, cdouble hz2i,
				 cdouble dt, cdouble nu,
				 cdouble bcs[][6],
				 cdouble tol, cuint max_iteration );


// sparse viscosity matrix
void viscosity_matrix_sparse( vector<tuple <uint, uint, double> >& L_sp,
							  vector<double>& val,
							  vector<uint>& col_ind,
							  vector<uint>& row_ptr,
							  double* F,
							  cuint nx, cuint ny, cuint nz,
							  cdouble hx, cdouble hy, cdouble hz,
							  cdouble hx2i, cdouble hy2i, cdouble hz2i,
							  cdouble dt, cdouble nu,
							  cdouble* u_bc,
							  cuint dir // direction of flow: u, v, or w?
							  );

// set load vector for implicit viscous solve
void viscosity_load_vector( double* F,  boost::multi_array<double, 3>& U);

#endif //VISCOSITY_H
