#include "jacobi.h"
#include "assemble.h"
#include "IO.h"

double convergence_check ( double** M,
						   double* U,
						   double* F,
						   double* R,
						   cuint n_dof)
{
	double E=0;
#pragma omp parallel for shared(R,M,U,F,E)
	for(int i=0; i<n_dof; i++){
		R[i] = 0.0;
		for(int j=0; j<n_dof; j++){
			R[i] -= M[i][j]*U[j];
		}
		R[i] += F[i];
		E += R[i]*R[i];
		
	}
	
	return E; 
}

// jacobi method
void jacobi( cdouble tol,
			 cuint max_iteration,
			 cuint n_dof,
			 double* u_new,
			 double* u_old,
			 double** M,
			 double* F,
			 double& Er,
			 double* R)
{
	// iteration counter
	int ct = 0;

	while(Er>tol && ct<max_iteration){
		for(int i=0;i<n_dof;i++)
			u_old[i]=u_new[i];
		
#pragma omp parallel for shared(M,F,u_old,u_new)
		for(int i=0; i<n_dof; i++){
			double S=0;
			for(int j=0; j<n_dof; j++){
				if(i!=j)
					S += M[i][j]*u_old[j]; 
			}
			// if(M[i][i]==0) cout<<"zero "<<i<<endl;
			// if(F[i] != F[i]) cout<<F[i]<<" "<<i<<endl;
			u_new[i] = 1/M[i][i] * (F[i] - S);
			// cout<<u_new[i]<<endl;
		}

		Er = convergence_check(M, u_new, F, R, n_dof);
		// cout<<"Er: "<<Er<<endl;
		ct++;
	}

	if(max_iteration==0) 		
		Er = convergence_check(M, u_new, F, R, n_dof);
	
	
	return;
}

// sparse jacobi method
/*
void jacobi_sparse( cdouble tol,
					cuint max_iteration,
					cuint n_dof,
					double* U,
					double* U_tmp,
					double** M,
					double* F,
					double& Er,
					double* R)
{
	// iteration counter
	int ct = 0;

	while(Er>tol && ct<max_iteration){
		for(int i=0;i<n_dof;i++)
			U_tmp[i]=U[i];
		
#pragma omp parallel for shared(M,F,U_tmp,U)
		while(n<M.size()){
			double S=get<2>(M[n]);
			cuint i= get<0>(M[n]);
			cuint j= get<1>(M[n]);
			while(get<0>(M[n])==i){

			}
		}
		
		
		for(int i=0; i<n_dof; i++){
			double S=0;
			for(int j=0; j<n_dof; j++){
				if(i!=j)
					S += M[i][j]*U_tmp[j]; 
			}
			// if(M[i][i]==0) cout<<"zero "<<i<<endl;
			// if(F[i] != F[i]) cout<<F[i]<<" "<<i<<endl;
			U[i] = 1/M[i][i] * (F[i] - S);
			// cout<<U[i]<<endl;
		}

		Er = convergence_check(M, U, F, R, n_dof);
		// cout<<"Er: "<<Er<<endl;
		ct++;
	}

	if(max_iteration==0) 		
		Er = convergence_check(M, U, F, R, n_dof);
	
	
	return;
}
*/

// multigrid v-cycle
double* v_cycle( uint n_dof, cuint I, cuint J, cuint K,
				 cdouble dx2i, cdouble dy2i, cdouble dz2i,
				 cdouble tol, cuint max_iteration, cuint pre_smooth_iteration,
				 cdouble width, cdouble length, cdouble height,
				 cuint level, cuint max_level,
				 double* F,
				 double& Er)
{
	cout<<"level: "<<level<<" n_dof: "<<n_dof<<endl;
	// for global constraint
	n_dof = n_dof+1;
	// initialize finite difference matrix (+1 for global constraint)
	double** M = new double*[n_dof];
	for(int n = 0; n < (n_dof); n++)
		M[n] = new double[n_dof];
	// initialize 
#pragma omp parallel for shared(n_dof, M)
	for(int i=0; i<n_dof; i++)
		for(int j=0; j<n_dof; j++)
			M[i][j] = 0;
	
	// create finite difference matrix
	cout<<"create finite difference matrix"<<endl;
	fd_matrix(M, I,J,K, dx2i, dy2i, dz2i, n_dof);

	// construct load vector
	// load vector is created only at the level 0
	if(level==0){
		F = new double[n_dof];
		cout<<"create load vector"<<endl;
		load_vector(F, n_dof, I,J,K );
	}

	// set dirichlet boundary conditions
	// unsigned int n_bd=boundary_conditins(n_dof, I, J, K, M, F);	
	// cout<<"number of boundary nodes = "<<n_bd<<endl;

	// cout<<"save matrix and vector"<<endl;
	// char matrix_file[100];
	// char vector_file[100];
	// sprintf(matrix_file, "matrix_%i.dat", level);
	// sprintf(vector_file, "vector_%i.dat", level);
	// if(write_matrix(n_dof,n_dof,M,matrix_file))
	// 	cout<<"write_matrix fail"<<endl;
	// if(write_vector(n_dof,F,vector_file)) cout<<"write_vector fail"<<endl;
	
	// construct solution vector
	double* U = new double[n_dof];
	double* U_tmp = new double[n_dof];
	// initial guess
#pragma omp parallel for shared(U, U_tmp)
	for(int n=0; n<n_dof; n++){
	    U[n] = 0.0;
	    U_tmp[n] = 0.0;
    }

	// residual and error
	double* R = new double[n_dof];

	// perform pre-smoothing and compute residual
	cout<<"pre-smoothing "<<pre_smooth_iteration<<" times"<<endl;
	jacobi(tol, pre_smooth_iteration, n_dof, U, U_tmp, M, F, Er, R);
		
	// restriction of residual on coarse grid
	double* F_coar;
		
	// Restrict the residual
	cuint I_coar = (I)/2;
	cuint J_coar = (J)/2;
	cuint K_coar = (K)/2;
	uint n_dof_coar = I_coar*J_coar*K_coar; // +1 for global constraint
	n_dof_coar +=1;
	F_coar = new double[n_dof_coar];

	// mesh size (-1) ignored because of periodic domain
	cdouble dx_coar = width/(I_coar);
	cdouble dy_coar = length/(J_coar);
	cdouble dz_coar = height/(K_coar);
	
	// inverse of square of mesh sizes
	cdouble dx2i_coar = 1.0/(dx_coar*dx_coar);
	cdouble dy2i_coar = 1.0/(dy_coar*dy_coar);
	cdouble dz2i_coar = 1.0/(dz_coar*dz_coar);
		
	// restric residual to the coarrse grid
	cout<<"restriction"<<endl;
	restriction( R, F_coar, I, J, K, I_coar, J_coar, K_coar);
	F_coar[n_dof_coar-1] = -R[n_dof-1]/n_dof*n_dof_coar; // global constraint
	
	// construct solution vector on coarse grid
	double* U_coar = new double[n_dof_coar];
	double* U_coar_tmp = new double[n_dof_coar];
	
	// if the grid is coarsest
	if( level==max_level){
		cout<<"level: "<<level+1<<" n_dof: "<<n_dof_coar<<endl;

		// initial guess
#pragma omp parallel for shared(U_coar, U_coar_tmp)
		for(int n=0; n<n_dof_coar; n++){
			U_coar[n] = 0.0;
			U_coar_tmp[n] = 0.0;
		}

		// initialize finite difference matrix (+1 for global constraint)
		double** M_coar = new double*[n_dof_coar];
		for(int n = 0; n < (n_dof_coar); n++)
			M_coar[n] = new double[n_dof_coar];
		// initialize 
#pragma omp parallel for shared(M_coar)
		for(int i=0; i<n_dof_coar; i++)
			for(int j=0; j<n_dof_coar; j++)
				M_coar[i][j] = 0;
		// create finite difference matrix
		cout<<"create finite difference matrix"<<endl;
		fd_matrix(M_coar, I_coar,J_coar,K_coar,
				  dx2i_coar, dy2i_coar, dz2i_coar, n_dof_coar);

		// residual on coarse grid
		double* R_coar = new double[n_dof_coar];
		
		// exact Jacobi method
		jacobi(tol, max_iteration, n_dof_coar, U_coar, U_coar_tmp,
			   M_coar, F_coar, Er, R_coar);

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
		U_coar = v_cycle( n_dof_coar-1, I_coar, J_coar, K_coar,
						  dx2i_coar, dy2i_coar, dz2i_coar,
						  tol, max_iteration, pre_smooth_iteration,
						  width, length, height, level+1, max_level, F_coar, Er );
		
		cdouble dx_coar = width/(I_coar);
		cdouble dy_coar = length/(J_coar);
		cdouble dz_coar = height/(K_coar);

	  
		// write_results( U_coar,
		// 			   n_dof_coar,
		// 			   I_coar, J_coar, K_coar,
		// 			   dx_coar, dy_coar, dz_coar, level);
		 
	}

	// fine grid (-1 ignored due to periodic domain)
	// cuint I_fine = I*2;
	// cuint J_fine = J*2;
	// cuint K_fine = K*2;
	// cuint n_dof_fine = I_fine*J_fine*K_fine+1;
	double* E = new double[n_dof];
	interpolation(U_coar, E, I_coar,J_coar,K_coar, I, J, K);

	// correct the fine grid approximation
#pragma omp parallel for shared(U,E)
	for(int i=0; i<n_dof; i++){
		// cout<<i<<" "<<U[i]<<" "<<E[i]<<" "<<E[i]/U[i]<<endl;
		U[i] += E[i];
	}

	// perform post-smoothing and compute residual
	uint post_smooth_iteration;
	if(level==0)
		post_smooth_iteration=max_iteration;
	else
		post_smooth_iteration=pre_smooth_iteration;

	cout<<"post-smoothing "<<post_smooth_iteration<<" times on level "
		<<level<<endl;
	jacobi(tol, post_smooth_iteration, n_dof, U, U_tmp, M, F, Er, R);

	// cleanup
	for(int n = 0; n< n_dof; n++) {
		delete[] M[n];
	}
	delete[] M;

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
#pragma omp parallel for shared(R, R_new) private(nei)
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

#pragma omp parallel for shared(U, U_fine) private(box_old, box_fine)
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
double* v_cycle_0( uint n_dof, cuint I, cuint J, cuint K,
				 cdouble dx2i, cdouble dy2i, cdouble dz2i,
				 cdouble tol, cuint max_iteration, cuint pre_smooth_iteration,
				 cdouble width, cdouble length, cdouble height,
				 cuint level, cuint max_level,
				 double* F,
				 double& Er)
{
	cout<<"level: "<<level<<" n_dof: "<<n_dof<<endl;
	// for global constraint
	n_dof = n_dof+1;
	// initialize finite difference matrix (+1 for global constraint)
	double** M = new double*[n_dof];
	for(int n = 0; n < (n_dof); n++)
		M[n] = new double[n_dof];
	// initialize 
#pragma omp parallel for shared(n_dof, M)
	for(int i=0; i<n_dof; i++)
		for(int j=0; j<n_dof; j++)
			M[i][j] = 0;
	
	// create finite difference matrix
	cout<<"create finite difference matrix"<<endl;
	fd_matrix(M, I,J,K, dx2i, dy2i, dz2i, n_dof);

	// construct load vector
	// load vector is created only at the level 0
	if(level==0){
		F = new double[n_dof];
		cout<<"create load vector"<<endl;
		load_vector(F, n_dof, I,J,K );
	}

	// set dirichlet boundary conditions
	// unsigned int n_bd=boundary_conditins(n_dof, I, J, K, M, F);	
	// cout<<"number of boundary nodes = "<<n_bd<<endl;

	// cout<<"save matrix and vector"<<endl;
	// char matrix_file[100];
	// char vector_file[100];
	// sprintf(matrix_file, "matrix_%i.dat", level);
	// sprintf(vector_file, "vector_%i.dat", level);
	// if(write_matrix(n_dof,n_dof,M,matrix_file))
	// 	cout<<"write_matrix fail"<<endl;
	// if(write_vector(n_dof,F,vector_file)) cout<<"write_vector fail"<<endl;
	
	// construct solution vector
	double* U = new double[n_dof];
	double* U_tmp = new double[n_dof];
	// initial guess
	for(int n=0; n<n_dof; n++){
	    U[n] = 0.0;
	    U_tmp[n] = 0.0;
    }

	// residual and error
	double* R = new double[n_dof];

	// perform full jacobi iteration
	cout<<"jacobi method "<<max_iteration<<" times"<<endl;
	jacobi(tol, max_iteration, n_dof, U, U_tmp, M, F, Er, R);
		
	// cleanup
	for(int n = 0; n< n_dof; n++) {
		delete[] M[n];
	}
	delete[] M;

	if (level==0)
		delete[] F;

	delete[] U_tmp;
	delete[] R;
	
	return U;
}
