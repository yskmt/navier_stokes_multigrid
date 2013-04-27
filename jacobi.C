#include "jacobi.h"
#include "assemble.h"

double convergence_check ( double** M,
						   double* U,
						   double* F,
						   double* R,
						   const int n_dof)
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
			 double& E,
			 double* R)
{
	// iteration counter
	int ct = 0;
	while(E>tol && ct<max_iteration){
		for(int i=0;i<n_dof;i++)
			u_old[i]=u_new[i];
		cout<<E<<endl;
		
#pragma omp parallel for shared(M,F,u_old,u_new)
		for(int i=0; i<n_dof; i++){
			double S=0;
			for(int j=0; j<n_dof; j++){
				if(i!=j)
					S += M[i][j]*u_old[j]; 
			}
			
			u_new[i] = 1/M[i][i] * (F[i] - S);
			// cout<<u_new[i]<<endl;
		}
		// check convergence
		// cout<<"conv check"<<endl;
		// if(!(ct%10))
		E = convergence_check(M, u_new, F, R, n_dof);
		ct++;
	}
	
	return;
}

void v_cycle( cuint n_dof, cuint I, cuint J, cuint K,
			  cdouble dx2i, cdouble dy2i, cdouble dz2i,
			  cdouble tol, cuint max_iteration,
			  cdouble width, cdouble length, cdouble height,
			  cuint level, cuint max_level,
			  double* F )
{		
	// initialize finite difference matrix
	double** M = new double*[n_dof];
	for(int n = 0; n < (n_dof); n++)
		M[n] = new double[n_dof];

	// create finite difference matrix
	fd_matrix(M, I,J,K, dx2i, dy2i, dz2i);
	
	// construct load vector
	// double* F;
	// load vector is created only at the level 0
	if(level==0){
		F = new double[n_dof];
		load_vector(F, n_dof, I,J,K);
	}
	// else the load vector is taken from residual
	// else{
	// 	F = R;
	// }

	// set boundary conditions
	unsigned int n_bd=boundary_conditins(n_dof, I, J, K, M);
	
	// cout<<"number of boundary nodes = "<<n_bd<<endl;

	// save matrix and vector
	// if(write_matrix(n_dof,n_dof,M)) cout<<"write_matrix fail"<<endl;
	// if(write_vector(n_dof,F)) cout<<"write_vector fail"<<endl;
	
	// construct solution vector
	double* u_new = new double[n_dof];
	double* u_old = new double[n_dof];
	// initial guess
	for(int n=0; n<n_dof; n++){
	    u_new[n] = 1.0;
	    u_old[n] = 1.0;
    }

	// residual and error
	double* R = new double[n_dof];
	double E=100000000;
	double* R_new;
	
	// Jacobi method
	cout<<"level: "<<level<<endl;
	jacobi(tol, max_iteration, n_dof, u_new, u_old, M, F, E, R);
	cout<<E<<endl;
	
	if( level<max_level){

		// Restrict the residual
		cuint I_new = (I-1)/2;
		cuint J_new = (J-1)/2;
		cuint K_new = (K-1)/2;
		cuint n_dof_new = I_new*J_new*K_new;
		R_new = new double[n_dof_new];

		// mesh size
		cdouble dx_new = width/(I_new-1);
		cdouble dy_new = length/(J_new-1);
		cdouble dz_new = height/(K_new-1);
	
		// inverse of square of mesh sizes
		cdouble dx2i_new = 1.0/(dx_new*dx_new);
		cdouble dy2i_new = 1.0/(dy_new*dy_new);
		cdouble dz2i_new = 1.0/(dz_new*dz_new);
	
		restriction( R, R_new, I, J, K);

		// for(int i=0; i<n_dof_new; i++)
		// 	cout<<R_new[i]<<endl;

		// v_cycle on the coarse grid
		v_cycle( n_dof_new, I_new, J_new, K_new,
				 dx2i_new, dy2i_new, dz2i_new,
				 tol, max_iteration,
				 width, length, height, level+1, max_level, R_new );

	}

	// cleanup
	for(int n = 0; n< n_dof; n++) {
		delete[] M[n];
	}
	delete[] M;

	if (level==0)
		delete[] F;

	delete[] u_new, u_old;
	delete[] R, R_new;
	
	return;
}
