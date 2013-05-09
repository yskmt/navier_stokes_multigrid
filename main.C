#include "jacobi.h"
#include "assemble.h"
#include "utils.h"
#include "IO.h"
#include "v_cycle.h"
#include "advection.h"
#include "viscosity.h"

// number of threads
uint nt;

int main( int argc, char** argv )
{
	// initialize constants
	cdouble nu = 100; // kinetic viscosity (mu/rho)
	double dt = 1e-2; //time step
	cdouble tf = 4.0; // final time
	
	// domain size
	cdouble lx = 1.0; 
	cdouble ly = 1.0;
	cdouble lz = 1.0;

	// domain cornders
	cdouble x_min = 0.0;
	cdouble y_min = 0.0;
	cdouble z_min = 0.0;
	cdouble x_max = x_min+lx;
	cdouble y_max = y_min+ly;
	cdouble z_max = z_min+lz;
	
	// number of gridpointts in each dimension
	nt=1;
	uint nx=16;
	uint ny=16;
	uint nz=16; // problem size (n_dof=n_size^3)
	uint max_level=0; // maximum v-cycle level
	if(argc>5){
		nt = atoi(argv[1]);
		max_level = atoi(argv[2]);
		nx = atoi(argv[3]);
		ny = atoi(argv[4]);
		nz = atoi(argv[5]);
	}
	else{
		cout<<"multigrid [# of threads] [max level] [I_size] [J_size] [K_size]"<<endl;
		return 0;
	}
	
	// number of nodes in each dimension
	// minimum size =3*3*3, should be 2^n+1: n=max_level-1
	// should be 2^n due to periodic domain

	cuint n_dof = nx*ny*nz;

	// number of time steps
	cuint nt = floor(tf/dt);
	// corrected time step size
	dt = tf/nt;

	// mesh size
	cdouble hx = lx/(nx);
	cdouble hy = ly/(ny);
	cdouble hz = lz/(nz);

	// inverse of square of mesh sizes
	cdouble hx2i = 1.0/(hx*hx);
	cdouble hy2i = 1.0/(hy*hy);
	cdouble hz2i = 1.0/(hz*hz);

	// set up initial conditions
	// interior values
	boost::multi_array<double, 3> U(boost::extents[nx-1][ny][nz]);
	boost::multi_array<double, 3> V(boost::extents[nx][ny-1][nz]);
	boost::multi_array<double, 3> W(boost::extents[nx][ny][nz-1]);
	
	// treat nonlinear (advection) terms
	advection(U,V,W, nx,ny,nz, hx, hy, hz, dt);

	cdouble bcs[3][6] = {{0,0,0,0,0,1.0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}};
	
	// implicitly solve viscosity terms
	viscosity( U, V, W, nx, ny, nz, hx, hy, hz,
			   hx2i, hy2i, hz2i,
			   dt, nu, bcs );
	
	// implicit viscosity terms
	// viscosity();
	
	/*
	// for jacobi method
	cdouble tol = 0.0001;
	cuint max_iteration = 10000;
	cuint pre_smooth_iteration = 10;

	double* F;
	double Er = tol*10;
	double* U;
	double start, end;

	start=omp_get_wtime();

// 	cout<<"fd_matrix_sparse"<<endl;
// 	vector<tuple <uint, uint, double> > M_sp;
// 	vector<double> val(M_sp.size(),0.0);
// 	vector<uint> col_ind(M_sp.size(), 0);
// 	vector<uint> row_ptr(1,0);
	
// 	fd_matrix_sparse(M_sp, val, col_ind, row_ptr,
// 					 I,J,K, hx2i, hy2i, hz2i, n_dof+1);
// 	cout<<"done"<<endl;

// 	double** M = new double*[n_dof+1];
// 	for(int n = 0; n < (n_dof+1); n++)
// 		M[n] = new double[n_dof+1];
// 	// initialize 
// #pragma omp parallel for shared(M)
// 	for(int i=0; i<n_dof+1; i++)
// 		for(int j=0; j<n_dof+1; j++)
// 			M[i][j] = 0;	
// 	fd_matrix(M, I,J,K, hx2i, hy2i, hz2i, n_dof+1);
// 	char file_name[]="test_matrix.dat";
// 	if(write_matrix(n_dof+1,n_dof+1,M, file_name))

	
	if(max_level==0){
		U = v_cycle_0( n_dof, nx, ny, nz,
				   hx2i, hy2i,  hz2i,
				   tol, max_iteration, pre_smooth_iteration,
					   lx, lz, lz, 0, max_level-1, F, Er );
	}
	else{
		U = v_cycle( n_dof, nx, ny, nz,
							 hx2i, hy2i,  hz2i,
							 tol, max_iteration, pre_smooth_iteration,
					 lx, ly, lz, 0, max_level, F, Er );
	}
	end=omp_get_wtime();
	cout<<"wall clock time = "<<end-start<<endl;
	cout<<"final error: "<<sqrt(Er)<<endl;

	
	// for(int i=0; i<n_dof; i++)
	// 	cout<<u_new[i]<<endl;
			
	write_results( U,
				   n_dof,
				   nx, ny, nz, hx, hy, hz, 100);
	
	// delete[] U;

	*/
	
	return 0;
}
