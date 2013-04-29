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
			R[i] += M[i][j]*U[j];
		}
		R[i] -= F[i];
		E += R[i]*R[i];
		
	}
	
	return E; 
}

// jacobi method
void jacobi( cdouble tol, const int max_iteration,
			 const unsigned int n_dof,
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
			// if(F[i] != F[i]) cout<<F[i]<<" "<<i<<endl;
			u_new[i] = 1/M[i][i] * (F[i] - S);
			// cout<<u_new[i]<<endl;
		}
		// check convrgence
		// cout<<"conv check"<<endl;
		// if(!(ct%10))

		Er = convergence_check(M, u_new, F, R, n_dof);
		cout<<"Er: "<<Er<<endl;
		ct++;
	}
	
	return;
}

// multigrid v-cycle
double* v_cycle( uint n_dof, cuint I, cuint J, cuint K,
				 cdouble dx2i, cdouble dy2i, cdouble dz2i,
				 cdouble tol, cuint max_iteration,
				 cdouble width, cdouble length, cdouble height,
				 cuint level, cuint max_level,
				 double* F )
{
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
	cout<<"cireate finite difference matrix"<<endl;
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

	cout<<"save matrix and vector"<<endl;
	char matrix_file[100];
	char vector_file[100];
	sprintf(matrix_file, "matrix_%i.dat", level);
	sprintf(vector_file, "vector_%i.dat", level);
	if(write_matrix(n_dof,n_dof,M,matrix_file))
		cout<<"write_matrix fail"<<endl;
	if(write_vector(n_dof,F,vector_file)) cout<<"write_vector fail"<<endl;
	
	// construct solution vector
	double* U = new double[n_dof];
	double* u_old = new double[n_dof];
	// initial guess
	for(int n=0; n<n_dof; n++){
	    U[n] = 0.0;
	    u_old[n] = 0.0;
    }

	// residual and error
	double* R = new double[n_dof];
	double Er=tol*100;
	double* R_coa;

	if( level==max_level){
		// exact Jacobi method
		cout<<"level: "<<level<<" n_dof "<<n_dof<<endl;

		jacobi(tol, max_iteration, n_dof, U, u_old, M, F, Er, R);

		// cout<<"R"<<endl;
		// for(int i=0; i<n_dof; i++)
		// 	cout<<R[i]<<endl;
		
	}
	else{
		// inexact Jacobi method
		cout<<"level: "<<level<<" n_dof "<<n_dof<<endl;
		jacobi(tol, 5, n_dof, U, u_old, M, F, Er, R);
		
		// Restrict the residual
		cuint I_coa = (I+1)/2;
		cuint J_coa = (J+1)/2;
		cuint K_coa = (K+1)/2;
		cuint n_dof_coa = I_coa*J_coa*K_coa;
		R_coa = new double[n_dof_coa];

		// mesh size (-1) ignored because of periodic domain
		cdouble dx_coa = width/(I_coa);
		cdouble dy_coa = length/(J_coa);
		cdouble dz_coa = height/(K_coa);
	
		// inverse of square of mesh sizes
		cdouble dx2i_coa = 1.0/(dx_coa*dx_coa);
		cdouble dy2i_coa = 1.0/(dy_coa*dy_coa);
		cdouble dz2i_coa = 1.0/(dz_coa*dz_coa);

		// for(int i=0; i<n_dof; i++)
		// 	if(R[i]!=R[i]) cout<<i<<" is nan"<<endl;
		
		// restric residual to the coarse grid
		cout<<"restriction"<<endl;
		restriction( R, R_coa, I, J, K, I_coa, J_coa, K_coa);
		cout<<"done"<<endl;	
		
		// v_cycle on the coarse grid
		v_cycle( n_dof_coa, I_coa, J_coa, K_coa,
				 dx2i_coa, dy2i_coa, dz2i_coa,
				 tol, max_iteration,
				 width, length, height, level+1, max_level, R_coa );

	}

	cuint I_fine = I*2-1;
	cuint J_fine = J*2-1;
	cuint K_fine = K*2-1;
	cuint n_dof_fine = I_fine*J_fine*K_fine+1;
	double* E = new double[n_dof_fine];
	interpolation(U, E, I,J,K, I_fine, J_fine, K_fine);
	
	// cleanup
	for(int n = 0; n< n_dof; n++) {
		delete[] M[n];
	}
	delete[] M;

	if (level==0)
		delete[] F;

	delete[] u_old;
	delete[] R, R_coa;
	
	return U;
}

// 3D full weight restriction
void restriction( double* R, double* R_new, cuint I, cuint J, cuint K,
				  cuint I_new, cuint J_new, cuint K_new )
{	
	unsigned int nei[3][3][3];
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
				 uint nei[][3][3], cuint i, cuint j, cuint k, cuint t_new )
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
