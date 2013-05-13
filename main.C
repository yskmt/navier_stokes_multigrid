#include "jacobi.h"
#include "assemble.h"
#include "utils.h"
#include "IO.h"
#include "v_cycle.h"
#include "advection.h"
#include "viscosity.h"
#include "pressure.h"

// number of threads
uint nt;
// set up initial conditions
void initial_conditions( double* U,
						 double* V,
						 double* W,
						 double* P,
						 cuint nx, cuint ny, cuint nz )
{
	uint t;
	// set up U
	for(int i=0; i<nx-1; i++){
		for(int j=0; j<ny; j++){
			for(int k=0; k<nz; k++){
				three_d_to_one_d(i,j,k, nx-1, ny, t);
				U[t] = 0.0;
			}
		}
	}

	// set up V
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny-1; j++){
			for(int k=0; k<nz; k++){
				three_d_to_one_d(i,j,k, nx, ny-1, t);
				V[t] = 0.0;
			}
		}
	}

	// set up W
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			for(int k=0; k<nz-1; k++){
				three_d_to_one_d(i,j,k, nx, ny, t);
				W[t] = 0.0;
			}
		}
	}	
	
	return;
}

// main function!
int main( int argc, char** argv )
{
	// initialize constants
	cdouble nu = 100; // kinetic viscosity (mu/rho)
	double dt = 1.0; //time step
	cdouble tf = 1.0; // final time
	
	// domain size
	cdouble lx = 1.0; 
	cdouble ly = 1.0;
	cdouble lz = 1.0;

	// domain cornders
	cdouble xmin = 0.0;
	cdouble ymin = 0.0;
	cdouble zmin = 0.0;
	cdouble xmax = xmin+lx;
	cdouble ymax = ymin+ly;
	cdouble zmax = zmin+lz;
	
	// number of gridpointts in each dimension
	nt=1;
	uint nx=10;
	uint ny=10;
	uint nz=10; // problem size (n_dof=n_size^3)
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
	cuint n_u_dof = (nx-1)*ny*nz;
	cuint n_v_dof = (nx)*(ny-1)*nz;
	cuint n_w_dof = (nx)*ny*(nz-1);
	
	// number of time steps
	cuint nts = floor(tf/dt);
	// corrected time step size
	dt = tf/nts;

	// mesh size
	cdouble hx = lx/(nx);
	cdouble hy = ly/(ny);
	cdouble hz = lz/(nz);

	// inverse of square of mesh sizes
	cdouble hx2i = 1.0/(hx*hx);
	cdouble hy2i = 1.0/(hy*hy);
	cdouble hz2i = 1.0/(hz*hz);

	// interior values
	double* U = new double[(nx-1)*(ny)*(nz)];
	double* V = new double[(nx)*(ny-1)*(nz)];
	double* W = new double[(nx)*(ny)*(nz-1)];
	double* P = new double[n_dof];

	// set up initial conditions
	cout<<"setting up initial conditions..."<<endl;
	initial_conditions(U, V, W, P, nx, ny, nz);
	
	// for jacobi method
	cdouble tol = 0.01;
	cuint max_iteration = 10000000;
	cuint pre_smooth_iteration = 10;
	
	// boundary conditions
	// x0 xl y0 yl z0 zl
	cdouble bcs[3][6] = { {0,0,0,0,0,1}, {0,0,0,0,0,1}, {0,0,0,0,0,0}};
	// cdouble bcs[3][6] = { {1,1,1,1,1,1}, {0,0,0,0,0,0}, {0,0,0,0,0,0}};

	double start=omp_get_wtime();	
	for(int ts=0; ts<nts; ts++){
		cout<<"loop :"<<ts<<endl;
		
		// treat nonlinear (advection) terms
		cout<<"calculating advection terms..."<<endl;
		advection(U,V,W, nx,ny,nz, hx, hy, hz, dt, bcs);
	
		// implicitly solve viscosity terms
		double* Uss = new double[n_u_dof];
		double* Vss = new double[n_v_dof];
		double* Wss = new double[n_w_dof];
		cout<<"solving for viscosity terms..."<<endl;
		viscosity( U, V, W, Uss, Vss, Wss, nx, ny, nz, hx, hy, hz,
				   hx2i, hy2i, hz2i,
				   dt, nu, bcs,
				   tol, max_iteration );

		// cout<<"Uss"<<endl;
		// for(int i=0; i<n_u_dof; i++)
		// 	cout<<Uss[i]<<endl;
		// cout<<endl;

		// cout<<"Vss"<<endl;
		// for(int i=0; i<n_v_dof; i++)
		// 	cout<<Vss[i]<<endl;
		// cout<<endl;
		
		// cout<<"Wss"<<endl;
		// for(int i=0; i<n_w_dof; i++)
		// 	cout<<Wss[i]<<endl;
		// cout<<endl;
				
		// solve for pressure and update
		cout<<"solving for pressure..."<<endl;
		pressure( U,V,W, P, Uss, Vss, Wss, nx, ny, nz, bcs, hx, hy, hz,
				  hx2i, hy2i, hz2i, tol, max_iteration );


		// cout<<"U"<<endl;
		// for(int i=0; i<n_u_dof; i++){
		// 	cout<<U[i]<<endl;
		// }
		// cout<<"V"<<endl;
		// for(int i=0; i<n_v_dof; i++){
		// 	cout<<V[i]<<endl;
		// }
		// cout<<"W"<<endl;
		// for(int i=0; i<n_w_dof; i++){
		// 	cout<<W[i]<<endl;
		// }

		
		// write out the results
		cout<<"writing results..."<<endl;
		write_results( U, V, W, P, n_dof, nx, ny, nz,
					   xmin, ymin, zmin,
					   hx, hy, hz, ts, bcs);

	}
	double end=omp_get_wtime();
	cout<<"wall clock time: "<<end-start<<" with "<<nt<<" threads"<<endl;
	
	/*
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

