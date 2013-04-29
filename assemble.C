#include "assemble.h"
#include "utils.h"

// 2nd order stencil
void fd_matrix( double** M,
				cuint I, cuint J, cuint K,
				const double dx2i,
				const double dy2i,
				const double dz2i,
				cuint n_dof
				)
{
#pragma omp parallel for shared(M)
	for(int i=0; i<I; i++){
		for(int j=0; j<J; j++){
			for(int k=0; k<K; k++){
				unsigned int p,q;
				unsigned int t_011,t_111,t_211,t_101,t_121,t_110,t_112;
				three_d_to_one_d(i,j,k, I,J, t_111);
				if(i==0)
					three_d_to_one_d(I-1,j,k, I,J, t_011);
				else
					three_d_to_one_d(i-1,j,k, I,J, t_011);
				if(i==(I-1))
					three_d_to_one_d(0,j,k, I,J, t_211);
				else
					three_d_to_one_d(i+1,j,k, I,J, t_211);

				if(j==0)
					three_d_to_one_d(i,J-1,k, I,J, t_101);
				else
					three_d_to_one_d(i,j-1,k, I,J, t_101);
				if(j==(J-1))
					three_d_to_one_d(i,0,k, I,J, t_121);
				else
					three_d_to_one_d(i,j+1,k, I,J, t_121);
								
				if(k==0)
					three_d_to_one_d(i,j,K-1, I,J, t_110);
				else
					three_d_to_one_d(i,j,k-1, I,J, t_110);
				if(k==(K-1))
					three_d_to_one_d(i,j,0, I,J, t_112);
				else
					three_d_to_one_d(i,j,k+1, I,J, t_112);
				
				// I
				M[t_111][t_011] += dx2i;
				M[t_111][t_111] += -2*dx2i;
				M[t_111][t_211] += dx2i;
				// J
				M[t_111][t_101] += dy2i;
				M[t_111][t_111] += -2*dy2i;
				M[t_111][t_121] += dy2i;
				// K
				M[t_111][t_110] += dz2i;
				M[t_111][t_111] += -2*dz2i;
				M[t_111][t_112] += dz2i;

			}
		}
	}

	// cout<<"setting global constraint"<<endl;
	// global constraint
	for(int i=0; i<(n_dof); i++){
		M[i][n_dof-1] = 1;
		M[n_dof-1][i] = 1;
	}
	M[n_dof-1][n_dof-1] = n_dof;

}

void load_vector( double* F,
				  cuint n_dof,
				  cuint I,
				  cuint J,
				  cuint K
				  )
{
	// construct load vector
	#pragma omp parallel for shared(F)
	for(int n=0; n<n_dof-1; n++){
		unsigned int i,j,k;
		one_d_to_three_d( n, I, J, i, j, k);
	    F[n] = sin(double(i)/double(I)*2*pi); // solution=(2pi)^2*sin(2pi*x);

	    // F[n] = sin(double(i)/double(I)*2*pi) * sin(double(j)/double(J)*2*pi)
			// * sin(double(k)/double(K)*2*pi);
     }

	// global constraint
	F[n_dof-1] = 0;

}

// set dirichlet boudnary condition
// for periodic domain, it is sufficient to set only one point
int boundary_conditins( const unsigned int n_dof,
		const unsigned int I,
		const unsigned int J,
		const unsigned int K,
		double** M,
		double* F
		)
 {
	int n_bd=0;
	// boundary conditions
	uint t;
	three_d_to_one_d(0,0,0, I,J, t);		
	for(int n=0; n<n_dof; n++){
	    M[n][t]=0;
		M[t][n]=0;
	}
	M[t][t] = 1;
	F[t] = 1;
	
	// #pragma omp parallel for shared(M)
	// for(int i=0; i<I; i++){
	// 	for(int j=0; j<J; j++){
	// 		for(int k=0; k<K; k++){
	// 			if(i==0 || j==0 || k==0
	// 				|| i==(I-1) || j==(J-1) || k==(K-1) ){
	// 				n_bd++;
	// 				unsigned int t;
	// 				three_d_to_one_d(i,j,k, I,J, t);
					
	// 				for(int n=0; n<n_dof; n++){
	// 					M[n][t]=0;
	// 					M[t][n]=0;
	// 				}
	// 				M[t][t] = 1;
	// 			}
	// 		}
	// 	}
	// }

	return n_bd;
}
 
