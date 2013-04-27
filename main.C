#include "jacobi.h"
#include "assemble.h"
#include "utils.h"
#include "IO.h"

int main()
{
	// number of nodes in each dimension
	cuint I=7; // 2^n-1
	cuint J=7;
	cuint K=7;
	cuint n_dof = I*J*K;
	
	// domain size
	cdouble width = 1.0; 
	cdouble length = 1.0;
	cdouble height = 1.0;

	// domain cornders
	cdouble x_min = 0.0;
	cdouble y_min = 0.0;
	cdouble z_min = 0.0;
	cdouble x_max = x_min+width;
	cdouble y_max = y_min+length;
	cdouble z_max = z_min+height;

	// mesh size
	cdouble dx = width/(I-1);
	cdouble dy = length/(J-1);
	cdouble dz = height/(K-1);

	// inverse of square of mesh sizes
	cdouble dx2i = 1.0/(dx*dx);
	cdouble dy2i = 1.0/(dy*dy);
	cdouble dz2i = 1.0/(dz*dz);

	// for jacobi method
	cdouble tol = 0.01;
	cuint max_iteration = 100000;
	cuint max_level=0;

	cdouble start=omp_get_wtime();
	double* F;
	double* u = v_cycle( n_dof, I, J, K,
			 dx2i, dy2i,  dz2i,
			 tol, max_iteration,
			 width, length, height, 0, max_level, F );
	cdouble end=omp_get_wtime();

	cout<<"wall clock time = "<<end-start<<endl;
	
	// for(int i=0; i<n_dof; i++)
	// 	cout<<u_new[i]<<endl;

	write_results( u,
				   n_dof,
				   I, J, K, dx, dy, dz);
	
	delete[] u;
	
	return 0;
}
