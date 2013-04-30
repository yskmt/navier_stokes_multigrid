#include "jacobi.h"
#include "assemble.h"
#include "utils.h"
#include "IO.h"

int main()
{
	// number of nodes in each dimension
	// minimum size =3*3*3, should be 2^n+1: n=max_level-1
	// should be 2^n due to periodic domain
	cuint I=16;
	cuint J=16;
	cuint K=16;
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
	cdouble dx = width/(I);
	cdouble dy = length/(J);
	cdouble dz = height/(K);

	// inverse of square of mesh sizes
	cdouble dx2i = 1.0/(dx*dx);
	cdouble dy2i = 1.0/(dy*dy);
	cdouble dz2i = 1.0/(dz*dz);

	// for jacobi method
	cdouble tol = 0.01;
	cuint max_iteration = 10000;
	cuint pre_smooth_iteration = 10;
	cuint max_level=1;

	double* F;
	double Er = tol*10;
	double* U;
	double start, end;

	start=omp_get_wtime();
	cout<<"fd_matrix_sparse"<<endl;
	vector<tuple <uint, uint, double> > M_sp;
	vector<double> val(M_sp.size(),0.0);
	vector<uint> col_ind(M_sp.size(), 0);
	vector<uint> row_ptr(1,0);
	
	fd_matrix_sparse(M_sp, val, col_ind, row_ptr,
					 I,J,K, dx2i, dy2i, dz2i, n_dof+1);
	cout<<"done"<<endl;

	double** M = new double*[n_dof+1];
	for(int n = 0; n < (n_dof+1); n++)
		M[n] = new double[n_dof+1];
	// initialize 
#pragma omp parallel for shared(M)
	for(int i=0; i<n_dof+1; i++)
		for(int j=0; j<n_dof+1; j++)
			M[i][j] = 0;	
	fd_matrix(M, I,J,K, dx2i, dy2i, dz2i, n_dof+1);
	char file_name[]="test_matrix.dat";
	if(write_matrix(n_dof+1,n_dof+1,M, file_name))

	
	// if(max_level==0){
	// 	U = v_cycle_0( n_dof, I, J, K,
	// 			   dx2i, dy2i,  dz2i,
	// 			   tol, max_iteration, pre_smooth_iteration,
	// 			   width, length, height, 0, max_level-1, F, Er);
	// }
	// else{
	// 	U = v_cycle( n_dof, I, J, K,
	// 						 dx2i, dy2i,  dz2i,
	// 						 tol, max_iteration, pre_smooth_iteration,
	// 						 width, length, height, 0, max_level, F, Er );
	// }
	end=omp_get_wtime();
	cout<<"wall clock time = "<<end-start<<endl;
	cout<<"final error: "<<Er;

	
	// for(int i=0; i<n_dof; i++)
	// 	cout<<u_new[i]<<endl;
			
	// write_results( U,
	// 			   n_dof,
	// 			   I, J, K, dx, dy, dz, 100);
	
	// delete[] U;
	
	return 0;
}
