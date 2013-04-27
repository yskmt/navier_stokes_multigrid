#include "jacobi.h"
		
void three_d_to_one_d( const unsigned int i,
					  const unsigned int j,
					  const unsigned int k,
					  const unsigned int I,
					  const unsigned int J,
					  unsigned int& t )
{
	t=i + j*I + k*I*J;
}

void one_d_to_three_d( const unsigned int t,
					   const unsigned int I,
					   const unsigned int J,
					   unsigned int& i,
					   unsigned int& j,
					   unsigned int& k)
{
	k = t/(I*J);
	j = (t-k*I*J)/I;
	i = t-j*I - k*I*J;
}


// write out the sparse matrix
int write_matrix(const unsigned int P,
				 const unsigned int Q,
				 double** U)
{
	ofstream file_out;
	file_out.open ("matrix.dat");

	if(!file_out.is_open()){
		return 1;
	}

	for(int p=0; p<P; p++){
		for(int q=0; q<Q; q++){
			file_out<<U[p][q]<<" ";
		}
		file_out<<endl;
	}
		
	file_out.close();
	
	return 0;
}


// write out the sparse matrix
int write_vector( const unsigned int P,
				  double* F)
{
	ofstream file_out;
	file_out.open ("vector.dat");

	if(!file_out.is_open()){
		return 1;
	}

	for(int p=0; p<P; p++){
		file_out<<F[p]<<endl;
	}
		
	file_out.close();
	
	return 0;
}

// write out the results
int write_results( double* u,
				   const int n_dof,
				   const unsigned int I,
				   const unsigned int J,
				   const unsigned int K )
{
	// double*** results;    // 3D results definition;
	// // begin memory allocation
	// results = new double**[I];
	// for(int x = 0; x < I; ++x) {
	// 	results[x] = new double*[J];
	// 	for(int y = 0; y < J; ++y) {
	// 		results[x][y] = new double[K];
	// 		for(int z = 0; z < K; ++z) { // initialize the values to whatever you want the default to be
	// 			results[x][y][z] = 0;
	// 		}
	// 	}
	// }

	// write out the results now
	ofstream file_out;
	file_out.open ("results.vtk");
	if(!file_out.is_open()){
		return 1;
	}
	// file_out<<"x coord, y coord, z coord, scalar"<<endl;
	// unsigned int i,j,k;
	// for(int n=0; n<n_dof; n++){
	// 	one_d_to_three_d( n, I, J, i, j, k);
	// 	file_out<<j<<", "<<i<<", "<<k<<", "<<u[n]<<endl;
	// }

	
	// header
	file_out<<"# vtk DataFile Version 3.0"<<endl;
	file_out<<"3d poisson problem"<<endl
			<<"ASCII"<<endl
			<<"DATASET STRUCTURED_GRID"<<endl
			<<"DIMENSIONS "<<I<<" "<<J<<" "<<K<<endl
			<<"POINTS "<<n_dof<<" "<<"float"<<endl;
	
	unsigned int i,j,k;
	for(int n=0; n<n_dof; n++){
		one_d_to_three_d( n, I, J, i, j, k);
		file_out<<i<<" "<<j<<" "<<k<<endl;
	}
	
	file_out<<"POINT_DATA "<<n_dof<<endl;
	file_out<<"SCALARS u float 1"<<endl
			<<"LOOKUP_TABLE default"<<endl;
	for(int n=0; n<n_dof; n++){
		file_out<<u[n]<<endl;
	}
			
	
		
	file_out.close();

}

int main()
{
	// number of grids
	const unsigned int I=16;
	const unsigned int J=16;
	const unsigned int K=16;
	const unsigned int n_dof = I*J*K;
	
	const unsigned int P = sqrt(I*J*K);
	const unsigned int Q = P;

	// domain size
	const double width = 1.0; 
	const double length = 1.0;
	const double height = 1.0;

	// domain cornders
	const double x_min = 0.0;
	const double y_min = 0.0;
	const double z_min = 0.0;
	const double x_max = x_min+width;
	const double y_max = y_min+length;
	const double z_max = z_min+height;

	// mesh size
	const double dx = width/I;
	const double dy = length/J;
	const double dz = height/K;

	// inverse of square of mesh sizes
	const double dx2i = 1.0/(dx*dx);
	const double dy2i = 1.0/(dy*dy);
	const double dz2i = 1.0/(dz*dz);

	// for jacobi method
	const double tol = 0.01;
	const int max_iteration = 100000;
	
	// initialize finite difference matrix
	double** M = new double*[n_dof];
	for(int n = 0; n < (n_dof); n++)
		M[n] = new double[n_dof];

	const double start=omp_get_wtime();
	// create finite difference matrix
	#pragma omp parallel for shared(M)
	for(int i=1; i<I-1; i++){
		for(int j=1; j<J-1; j++){
			for(int k=1; k<K-1; k++){
				unsigned int p,q;
				unsigned int t_011,t_111,t_211,t_101,t_121,t_110,t_112;
				// I
				three_d_to_one_d(i-1,j,k, I,J, t_011);
				three_d_to_one_d(i,j,k, I,J, t_111);
				three_d_to_one_d(i+1,j,k, I,J, t_211);
				M[t_111][t_011] += dx2i;
				M[t_111][t_111] += -2*dx2i;
				M[t_111][t_211] += dx2i;

				// J
				three_d_to_one_d(i,j-1,k, I,J, t_101);
				// three_d_to_one_d(i,j,k, I,J, t_111);
				three_d_to_one_d(i,j+1,k, I,J, t_121);
				M[t_111][t_101] += dy2i;
				M[t_111][t_111] += -2*dy2i;
				M[t_111][t_121] += dy2i;

				// K
				three_d_to_one_d(i,j,k-1, I,J, t_110);
				// three_d_to_one_d(i,j,k, I,J, t_111);
				three_d_to_one_d(i,j,k+1, I,J, t_112);
				M[t_111][t_110] += dz2i;
				M[t_111][t_111] += -2*dz2i;
				M[t_111][t_112] += dz2i;
								
			}
		}
	}
	
	// construct load vector
	double* F = new double[n_dof];
	#pragma omp parallel for shared(F)
	for(int n=0; n<n_dof; n++){
		unsigned int i,j,k;
		one_d_to_three_d( n, I, J, i, j, k);
	    F[n] = sin(i/I*pi)*sin(j/J*pi)*sin(k/K*pi);
    }

	int n_bd=0;
	// boundary conditions
	#pragma omp parallel for shared(M)
	for(int i=0; i<I; i++){
		for(int j=0; j<J; j++){
			for(int k=0; k<K; k++){
				if(i==0 || j==0 || k==0
				   || i==(I-1) || j==(J-1) || k==(K-1) ){
					n_bd++;
					unsigned int t;
					three_d_to_one_d(i,j,k, I,J, t);

					for(int n=0; n<n_dof; n++){
						M[n][t]=0;
						M[t][n]=0;
					}
					M[t][t] = 1;
				}
				
			}
		}
	}
	cout<<"number of boundary nodes = "<<n_bd<<endl;
	// for(int p=0; p<P; p++){
	// 	for(int q=0; q<Q; q++){
	// 		cout<<M[p][q]<<" ";
	// 	}
	// 	cout<<endl;
	// }

	// save matrix and vector
	// if(write_matrix(n_dof,n_dof,M)) cout<<"write_matrix fail"<<endl;
	// if(write_vector(n_dof,F)) cout<<"write_vector fail"<<endl;

	
	// Jacobi method
	double E=100000000;

	// construct solution vector
	double* u_new = new double[n_dof];
	double* u_old = new double[n_dof];
	// initial guess
	for(int n=0; n<n_dof; n++){
	    u_new[n] = 1.0;
	    u_old[n] = 1.0;
    }

	jacobi(tol, max_iteration, n_dof, u_new, u_old, M, F, E);
	
	const double end=omp_get_wtime();

	cout<<"wall clock time = "<<end-start<<endl;
		
	// for(int i=0; i<n_dof; i++)
	// 	cout<<u_new[i]<<endl;

	write_results( u_new,
				   n_dof,
				   I,
				   J,
				   K );
	
	// cleanup
	for(int n = 0; n< n_dof; n++) {
		delete[] M[n];
	}
	delete[] M;
	delete[] F;
	
	return 0;
}
