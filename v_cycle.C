#include "v_cycle.h"
#include "pressure.h"

// multigrid v-cycle
double* v_cycle( uint n_dof, cuint nx, cuint ny, cuint nz,
				 cdouble hx, cdouble hy, cdouble hz,
				 cdouble hx2i, cdouble hy2i, cdouble hz2i,
				 cdouble tol, cuint max_iteration, cuint pre_smooth_iteration,
				 cdouble lx, cdouble ly, cdouble lz,
				 cuint level, cuint max_level,
				 double* F,
				 double& Er,
				 double* Uss, double* Vss, double* Wss,
				 cdouble bcs[][6]
				  )
{
	cout<<"level: "<<level<<" n_dof: "<<n_dof<<endl;

	// initialize finite difference matrix (+1 for global constraint)
// 	double** M = new double*[n_dof];
// 	for(int n = 0; n < (n_dof); n++)
// 		M[n] = new double[n_dof];
// 	// initialize 
// #pragma omp parallel for shared(n_dof, M)
// 	for(int i=0; i<n_dof; i++)
// 		for(int j=0; j<n_dof; j++)
// 			M[i][j] = 0;

	cout<<"fd_matrix_sparse"<<endl;
	vector<tuple <uint, uint, double> > M_sp;
	vector<double> val;
	vector<uint> col_ind;
	vector<uint> row_ptr(1,0);
	
	// create finite difference matrix
	cout<<"create finite difference matrix"<<endl;
	// build pressure matrix
	pressure_matrix( M_sp,
					 val, col_ind, row_ptr,
					 nx, ny, nz,
					 hx2i, hy2i, hz2i,
					 n_dof
					 );	
	
	// construct load vector
	// load vector is created only at the level 0
	if(level==0){
		F = new double[n_dof];
		cout<<"create load vector"<<endl;

		pressure_rhs(F, Uss, Vss, Wss, nx, ny, nz, bcs, hx, hy, hz);
		// load_vector(F, n_dof, I,J,K );
	}

	// cout<<"save matrix and vector"<<endl;
	// char matrix_file[100];
	char vector_file[100];
	// sprintf(matrix_file, "matrix_%i.dat", level);
	sprintf(vector_file, "vector_%i.dat", level);
	// if(write_matrix(n_dof,n_dof,M,matrix_file))
	// 	cout<<"write_matrix fail"<<endl;
	if(write_vector(n_dof,F,vector_file)) cout<<"write_vector fail"<<endl;
	
	// construct solution vector
	double* U = new double[n_dof];
	double* U_tmp = new double[n_dof];
	// initial guess
#pragma omp parallel for shared(U, U_tmp) num_threads(nt)
	for(int n=0; n<n_dof; n++){
	    U[n] = 0.0;
	    U_tmp[n] = 0.0;
    }

	// residual and error
	double* R = new double[n_dof];

	// perform pre-smoothing and compute residual
	cout<<"pre-smoothing "<<pre_smooth_iteration<<" times"<<endl;
	Er = tol*10;
	jacobi_sparse(tol, pre_smooth_iteration, n_dof, U, U_tmp,
				  val, col_ind, row_ptr, F, Er, R);
		
	// restriction of residual on coarse grid
	double* F_coar;
		
	// Restrict the residual
	cuint nx_coar = (nx)/2;
	cuint ny_coar = (ny)/2;
	cuint nz_coar = (nz)/2;
	uint n_dof_coar = nx_coar*ny_coar*nz_coar; 
	F_coar = new double[n_dof_coar];

	// mesh size 
	cdouble hx_coar = lx/(nx_coar);
	cdouble hy_coar = ly/(ny_coar);
	cdouble hz_coar = lz/(nz_coar);
	
	// inverse of square of mesh sizes
	cdouble hx2i_coar = 1.0/(hx_coar*hx_coar);
	cdouble hy2i_coar = 1.0/(hy_coar*hy_coar);
	cdouble hz2i_coar = 1.0/(hz_coar*hz_coar);
		
	// restric residual to the coarrse grid
	cout<<"restriction"<<endl;
	restriction( R, F_coar, nx, ny, nz, nx_coar, ny_coar, nz_coar);
	
	// construct solution vector on coarse grid
	double* U_coar = new double[n_dof_coar];
	double* U_coar_tmp = new double[n_dof_coar];
	
	// if the grid is coarsest
	if( level==max_level){
		cout<<"level: "<<level+1<<" n_dof: "<<n_dof_coar<<endl;

		// initial guess
#pragma omp parallel for shared(U_coar, U_coar_tmp) num_threads(nt)
		for(int n=0; n<n_dof_coar; n++){
			U_coar[n] = 0.0;
			U_coar_tmp[n] = 0.0;
		}

		vector<tuple <uint, uint, double> > M_sp_coar;
		vector<double> val_coar;
		vector<uint> col_ind_coar;
		vector<uint> row_ptr_coar(1,0);
		
		// create finite difference matrix
		cout<<"create finite difference matrix"<<endl;
		// fd_matrix_sparse(M_sp_coar, val_coar, col_ind_coar, row_ptr_coar,
		// 				 nx_coar,ny_coar,nz_coar,
		// 				 hx2i_coar, hy2i_coar, hz2i_coar, n_dof_coar );
		
		pressure_matrix( M_sp_coar, val_coar, col_ind_coar, row_ptr_coar,
						 nx_coar, ny_coar, nz_coar,
						 hx2i_coar, hy2i_coar, hz2i_coar,
						 n_dof_coar
						 );
		
		
		// residual on coarse grid
		double* R_coar = new double[n_dof_coar];
		
		// exact Jacobi method
		Er = tol*10;
		jacobi_sparse(tol, max_iteration, n_dof_coar, U_coar, U_coar_tmp,
					  val_coar, col_ind_coar, row_ptr_coar, F_coar,
					  Er, R_coar);
		
		// write_results( U_coar,
		// 			   n_dof_coar,
		// 			   I_coar, J_coar, K_coar,
		// 			   dx_coar, dy_coar, dz_coar, level);
		
		delete[] R_coar;
		
		// cout<<"R"<<endl;
		// for(int i=0; i<n_dof; i++)
		// 	cout<<R[i]<<endl;
		 
	}
	else{
		// v_cycle on the coarse grid
		U_coar = v_cycle( n_dof_coar, nx_coar, ny_coar, nz_coar,
						  hx_coar, hy_coar, hz_coar,
						  hx2i_coar, hy2i_coar, hz2i_coar,
						  tol, max_iteration, pre_smooth_iteration,
						  lx, ly, lz,
						  level+1, max_level,
						  F_coar, Er,
						  Uss, Vss, Wss,
						  bcs
						  );
		
		cdouble dx_coar = lx/(nx_coar);
		cdouble dy_coar = ly/(ny_coar);
		cdouble dz_coar = lz/(nz_coar);

		// // write partial results for test purpose
		// write_results( U_coar,
		// 			   n_dof_coar,
		// 			   I_coar, J_coar, K_coar,
		// 			   dx_coar, dy_coar, dz_coar, level);
		 
	}

	// interpolate to fine grid
	double* E = new double[n_dof];
	interpolation(U_coar, E, nx_coar,ny_coar,nz_coar, nx, ny, nz);

	// correct the fine grid approximation
#pragma omp parallel for shared(U,E) num_threads(nt)
	for(int i=0; i<n_dof; i++){
		// cout<<i<<" "<<U[i]<<" "<<E[i]<<" "<<E[i]/U[i]<<endl;
		U[i] += E[i];
	}

	// perform post-smoothing and compute residual
	uint post_smooth_iteration;
	// if(level==0)
		post_smooth_iteration=max_iteration;
	// else
		// post_smooth_iteration=( pre_smooth_iteration+1)*1000;

	cout<<"post-smoothing "<<post_smooth_iteration<<" times on level "
		<<level<<endl;
	// jacobi(tol, post_smooth_iteration, n_dof, U, U_tmp, M, F, Er, R);
	Er = tol*10;
	jacobi_sparse(tol, post_smooth_iteration, n_dof, U, U_tmp,
				  val, col_ind, row_ptr, F, Er, R);

	
	// cleanup
	if (level==0)
		delete[] F;

	delete[] U_tmp;
	delete[] R, F_coar;
	delete[] E;
	delete[] U_coar, U_coar_tmp;
	
	return U;
}

// 3D full weight restriction
void restriction( double* R, double* R_new, cuint I, cuint J, cuint K,
				  cuint I_new, cuint J_new, cuint K_new )
{	
	unsigned int nei[3][3][3];
#pragma omp parallel for shared(R, R_new) private(nei) num_threads(nt)
	for(int i=0; i<I; i+=2){
		for(int j=0; j<J; j+=2){
			for(int k=0; k<K; k+=2){
				get_neighbor( nei, i,j,k, I,J,K);
				// get new index
				uint t_new;
				three_d_to_one_d( i/2, j/2, k/2, I_new, J_new, t_new);
				coarse_map( R, R_new, nei, i,j,k, t_new);
			}
		}
	}
}

// map from fine to coarse grid
void coarse_map( double* R, double* R_new,
				 cuint nei[][3][3], cuint i, cuint j, cuint k, cuint t_new )
{
	// initialize to 0 (otherwise nan can )
	R_new[t_new] = 0.0;
	for(int p=0; p<3; p++){
		for(int q=0; q<3; q++){
			for(int r=0; r<3; r++){
				R_new[t_new] += R[nei[p][q][r]]*fw_stencil[p][q][r];
			}
		}
	}
}

// 3d trilinear interpolation
void interpolation( double* U, double* U_fine,
					cuint I, cuint J, cuint K,
					cuint I_fine, cuint J_fine, cuint K_fine)
{

	uint box_old[2][2][2];
	uint box_fine[2][2][2];

#pragma omp parallel for shared(U, U_fine) private(box_old, box_fine) num_threads(nt)
	for(int i=0; i<I; i++){
		for(int j=0; j<J; j++){
			for(int k=0; k<K; k++){
				// get the node nubmers of the old (coarse) box
				get_box(box_old, i,j,k, I, J, K);
				// get the node nubmbers of new (fine) box
				get_box(box_fine, i*2,j*2,k*2, I_fine, J_fine, K_fine);

				// map from coarse to fine grid
				fine_map(U, U_fine, box_old, box_fine);
			}
		}
	}

	// global constraints
	U_fine[I_fine*J_fine*K_fine] =
		U[I*J*K]/(I*J*K)*I_fine*J_fine*K_fine; // global constraint

	
	return;
}

void fine_map( double* U, double* U_new,
			   uint box_old[][2][2],
			   uint box_new[][2][2] )
{
	U_new[box_new[0][0][0]] = U[box_old[0][0][0]];
	U_new[box_new[1][0][0]] = (U[box_old[0][0][0]]
							   + U[box_old[1][0][0]])/2;
	U_new[box_new[0][1][0]] = (U[box_old[0][0][0]]
							   + U[box_old[0][1][0]])/2;
	U_new[box_new[0][0][1]] = (U[box_old[0][0][0]]
							   + U[box_old[0][0][1]])/2;
	U_new[box_new[1][1][0]] = (U[box_old[0][0][0]]
							   + U[box_old[1][0][0]]
							   + U[box_old[0][1][0]]
							   + U[box_old[1][1][0]])/4;
	U_new[box_new[1][0][1]] = (U[box_old[0][0][0]]
							   + U[box_old[1][0][0]]
							   + U[box_old[0][0][1]]
							   + U[box_old[1][0][1]])/4;
	U_new[box_new[0][1][1]] = (U[box_old[0][0][0]]
							   + U[box_old[0][1][0]]
							   + U[box_old[0][0][1]]
							   + U[box_old[0][1][1]])/4;
	U_new[box_new[1][1][1]] = ( U[box_old[0][0][0]]
								+ U[box_old[1][0][0]]
								+ U[box_old[0][1][0]]
								+ U[box_old[0][0][1]]
								+ U[box_old[0][1][1]]
								+ U[box_old[1][0][1]]
								+ U[box_old[1][1][0]]
								+ U[box_old[1][1][1]]
								)/8;

}


// 0_level v-cycle for testing purpose
double* v_cycle_0( double* Rp,
				   uint n_dof, cuint nx, cuint ny, cuint nz,
				   cdouble hx, cdouble hy, cdouble hz,
				   cdouble hx2i, cdouble hy2i, cdouble hz2i,
				   cdouble tol, cuint max_iteration, cuint pre_smooth_iteration,
				   cdouble width, cdouble length, cdouble height,
				   cuint level, cuint max_level,
				   double& Er,
				   double* Uss, double* Vss, double* Wss,
				   cdouble bcs[][6]				   
				   )
{
	cout<<"level: "<<level<<" n_dof: "<<n_dof<<endl;



	// load vector (extra +1 for global constraint) 
	double* Fp = new double[n_dof];

	// Lp
	vector<tuple <uint, uint, double> > Lp_sp;
	vector<double> Lp_val(Lp_sp.size(),0.0);
	vector<uint> Lp_col_ind(Lp_sp.size(), 0);
	vector<uint> Lp_row_ptr(1,0);		
	
	// build right hand side of pressure poisson equation
	pressure_rhs(Fp, Uss, Vss, Wss, nx, ny, nz, bcs, hx, hy, hz);

	// build pressure matrix
	pressure_matrix( Lp_sp,
					 Lp_val, Lp_col_ind, Lp_row_ptr,
					 nx, ny, nz,
					 hx2i, hy2i, hz2i,
					 n_dof
					 );	

	// solve dicrete poisson equation: Lp\Fp
	// construct solution vector
	double* P = new double[n_dof];
	double* P_tmp = new double[n_dof];
	// initial guess
#pragma omp parallel for shared(P, P_tmp) num_threads(nt)
	for(int n=0; n<n_dof; n++){
	    P[n] = 0.0;
	    P_tmp[n] = 0.0;
    }
	
	// jacobi iteration
	jacobi_sparse(tol, max_iteration, n_dof, P, P_tmp,
				  Lp_val, Lp_col_ind, Lp_row_ptr, Fp, Er, Rp);


	return P;
}
