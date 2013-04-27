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
				   cuint n_dof,
				   cuint I,
				   cuint J,
				   cuint K,
				   cdouble dx,
				   cdouble dy,
				   cdouble dz
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
	cuint I=15; // 2^n-1
	cuint J=15;
	cuint K=15;
	cuint n_dof = I*J*K;
	
	// domain size
	cdouble width = 1.0; 
	cdouble length = 1.0;
	cdouble height = 1.0;

	// domain cornders
	cdouble x_min = 0.0;
	cdouble y_min = 0.0;
	cdouble z_min = 0.0;
	cdouble x_max = x_min+width;
	cdouble y_max = y_min+length;
	cdouble z_max = z_min+height;

	// mesh size
	cdouble dx = width/(I-1);
	cdouble dy = length/(J-1);
	cdouble dz = height/(K-1);

	// inverse of square of mesh sizes
	cdouble dx2i = 1.0/(dx*dx);
	cdouble dy2i = 1.0/(dy*dy);
	cdouble dz2i = 1.0/(dz*dz);

	// for jacobi method
	cdouble tol = 0.01;
	cuint max_iteration = 100000;
	cuint max_level=2;

	cdouble start=omp_get_wtime();
	double* F;
	v_cycle( n_dof, I, J, K,
			 dx2i, dy2i,  dz2i,
			 tol, max_iteration,
			 width, length, height, 0, max_level, F );
	cdouble end=omp_get_wtime();

	cout<<"wall clock time = "<<end-start<<endl;
	
	// for(int i=0; i<n_dof; i++)
	// 	cout<<u_new[i]<<endl;

	// write_results( u_new,
	// 			   n_dof,
	// 			   I, J, K, dx, dy, dz);
	
	return 0;
}
