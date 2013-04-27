#include "jacobi.h"

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


 void jacobi( const double tol, const int max_iteration,
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
	cout<<"Jacobi"<<endl;
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
