#include "jacobi.h"
#include "assemble.h"
#include "utils.h"

// write out the sparse matrix
int write_matrix(cuint P,
				 cuint Q,
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
int write_vector( cuint P,
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
				   cuint I,
				   cuint J,
				   cuint K,
				   const double dx,
				   const double dy,
				   const double dz
				   )
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
		file_out<<i*dx<<" "<<j*dy<<" "<<k*dz<<endl;
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
	// number of nodes in each dimension
	cuint I=7; // 2^n-1
	cuint J=7;
	cuint K=7;
	cuint n_dof = I*J*K;
	
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
	const double dx = width/(I-1);
	const double dy = length/(J-1);
	const double dz = height/(K-1);

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
	fd_matrix(M, I,J,K, dx2i, dy2i, dz2i);
	
	// construct load vector
	double* F = new double[n_dof];
	load_vector(F, n_dof, I,J,K);

	// set boundary conditions
	unsigned int n_bd=boundary_conditins(n_dof, I, J, K, M);
	

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

	// residual
	double* R = new double[n_dof];

	jacobi(tol, max_iteration, n_dof, u_new, u_old, M, F, E, R);

	cuint I_new = (I-1)/2;
	cuint J_new = (J-1)/2;
	cuint K_new = (K-1)/2;
	cuint n_dof_new = I_new*J_new*K_new;
	double* R_new = new double[n_dof_new];
	restriction( R, R_new, I, J, K);

	// for(int i=0; i<n_dof_new; i++)
	// 	cout<<R_new[i]<<endl;
	
	const double end=omp_get_wtime();

	cout<<"wall clock time = "<<end-start<<endl;
		
	// for(int i=0; i<n_dof; i++)
	// 	cout<<u_new[i]<<endl;

	write_results( u_new,
				   n_dof,
				   I, J, K, dx, dy, dz);
	
	// cleanup
	for(int n = 0; n< n_dof; n++) {
		delete[] M[n];
	}
	delete[] M;
	delete[] F;
	
	return 0;
}
