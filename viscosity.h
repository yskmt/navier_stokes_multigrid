// implicit viscosity contributions
#ifndef VISCOSITY_H
#define VISCOSITY_H

#include "utils.h"
#include "assemble.h"
#include "msort.h"

using namespace std;

// implicit solve viscosity
// implicitly solve viscosity
void viscosity( cuint nx, cuint ny, cuint nz,
					cdouble hx2i, cdouble hy2i, cdouble hz2i,
					cdouble dt, cdouble nu );


// sparse viscosity matrix
void viscosity_matrix_sparse( vector<tuple <uint, uint, double> >& L_sp,
							  vector<double>& val,
							  vector<uint>& col_ind,
							  vector<uint>& row_ptr,
							  cuint nx, cuint ny, cuint nz,
							  cdouble hx2i,
							  cdouble hy2i,
							  cdouble hz2i,
							  cuint n_dof,
							  cuint dir
							  );

#endif //VISCOSITY_H
