#include "jacobi.h"
#include "assemble.h"
#include "IO.h"

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
	cdouble tol2 = tol*tol;
	
	while(Er>tol2 && ct<max_iteration){
		for(int i=0;i<n_dof;i++)
			u_old[i]=u_new[i];
		
#pragma omp parallel for shared(M,F,u_old,u_new) num_threads(nt)
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

	cout<<"convergence reached after "<<ct<<"iterations"<<endl;
	
	return;
}

// sparse jacobi method
void jacobi_sparse( cdouble tol,
					cuint max_iteration,
					cuint n_dof,
					double* U,
					double* U_tmp,
					const vector<double>& val,
					const vector<uint>& col_ind,
					const vector<uint>& row_ptr,
					double* F,
					double& Er,
					double* R)
{
	// iteration counter
	int ct = 0;
	cdouble tol2 = tol*tol;
	
	while(Er>tol2 && ct<max_iteration){
		for(int i=0;i<n_dof;i++)
			U_tmp[i]=U[i];
		
#pragma omp parallel for num_threads(nt) shared(F,U_tmp,U) 
 		for(int i=0; i<row_ptr.size()-1; i++){
			double S=0;
			double T=0;
			for(int j=row_ptr[i]; j<row_ptr[i+1]; j++){
				// cout<<"U_tmp "<<U_tmp[col_ind[j]]<<endl;
				if(i!=col_ind[j])
					S += val[j]*U_tmp[col_ind[j]];
				else{ // get diagonal element
					T = val[j];					
				}
			}
			U[i] = 1/T * (F[i]-S);
		}
		
		Er = convergence_check_sparse(val, col_ind, row_ptr, U, F, R, n_dof);
		// cout<<"i: "<<ct<<" Er: "<<Er<<endl;
		ct++;
	}

	if(max_iteration==0) 		
		// Er = convergence_check(M, U, F, R, n_dof);
		Er = convergence_check_sparse(val, col_ind, row_ptr, U, F, R, n_dof);

	
	cout<<"convergence reached after "<<ct<<" iterations"<<endl;
	
	
	return;
}

double convergence_check ( double** M,
						   double* U,
						   double* F,
						   double* R,
						   cuint n_dof
						   )
{
	double E=0;
#pragma omp parallel for shared(R,M,U,F,E) num_threads(nt)
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

double convergence_check_sparse ( const vector<double>& val,
								  const vector<uint>& col_ind,
								  const vector<uint>& row_ptr,
								  double* U,
								  double* F,
								  double* R,
								  cuint n_dof)
{
	double E=0;
#pragma omp parallel for shared(R,val,col_ind,row_ptr,U,F,E) num_threads(nt)
	for(int i=0; i<row_ptr.size()-1; i++){
		R[i] = 0.0;
		for(int j=row_ptr[i]; j<row_ptr[i+1]; j++){
			R[i] -= val[j]*U[col_ind[j]];
		}
		R[i] += F[i];
		E += R[i]*R[i];
	}
	
	return E; 
}
