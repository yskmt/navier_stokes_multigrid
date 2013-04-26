#include <iostream>
#include <cmath>
#include <fstream>
#include <omp.h>

using namespace std;

int three_d_to_one_d( const unsigned int i,
					  const unsigned int j,
					  const unsigned int k,
					  const unsigned int I,
					  const unsigned int J)
{
	return i + j*I + k*I*J;
}

void one_d_to_two_d( const unsigned int t,
					 const unsigned int P,
					 unsigned int& p,
					 unsigned int& q )
{
	p = t%P;
	q = t/P;
}

void three_d_to_one_d( const unsigned int i,
					   const unsigned int j,
					   const unsigned int k,
					   const unsigned int I,
					   const unsigned int J,
					   const unsigned int P,
					   unsigned int& p,
					   unsigned int& q)
{
	one_d_to_two_d(three_d_to_one_d(i,j,k, I, J), P, p, q);
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


int main()
{
	// number of grids
	const unsigned int I=128;
	const unsigned int J=128;
	const unsigned int K=128;

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

	cout<<dx2i<<endl;
	
	// initialize finite difference matrix
	double** U = new double*[P];
	for(int q = 0; q < P; ++q)
		U[q] = new double[Q];

	const double start=omp_get_wtime();
	// create finite difference matrix
	#pragma omp parallel for shared(U)
	for(int i=1; i<I-1; i++){
		for(int j=1; j<J-1; j++){
			for(int k=1; k<K-1; k++){
				unsigned int p,q;
				// I
				three_d_to_one_d(i-1,j,k, I,J,P, p,q);
				// if(p>=0 && p<P && q>=0 && q<Q)
				U[p][q] += dx2i;

				three_d_to_one_d(i,j,k, I,J,P, p,q);
				// if(p>=0 && p<P && q>=0 && q<Q)
				U[p][q] += -2*dx2i;
				
				three_d_to_one_d(i+1,j,k, I,J,P, p,q);
				// if(p>=0 && p<P && q>=0 && q<Q)
				U[p][q] += dx2i;

				// J
				three_d_to_one_d(i,j-1,k, I,J,P, p,q);
				// if(p>=0 && p<P && q>=0 && q<Q)
				U[p][q] += dy2i;
				
				three_d_to_one_d(i,j,k, I,J,P, p,q);
				// if(p>=0 && p<P && q>=0 && q<Q)
				U[p][q] += -2*dy2i;
				
				three_d_to_one_d(i,j+1,k, I,J,P, p,q);
				// if(p>=0 && p<P && q>=0 && q<Q)
				U[p][q] += dy2i;

				// K
				three_d_to_one_d(i,j,k-1, I,J,P, p,q);
				// if(p>=0 && p<P && q>=0 && q<Q)
				U[p][q] += dz2i;

				three_d_to_one_d(i,j,k, I,J,P, p,q);
				// if(p>=0 && p<P && q>=0 && q<Q)
				U[p][q] += -2*dz2i;
				
				three_d_to_one_d(i,j,k+1, I,J,P, p,q);
				// if(p>=0 && p<P && q>=0 && q<Q)
				U[p][q] += dz2i;
				
			}
		}
	}
	const double end=omp_get_wtime();

	cout<<end-start<<endl;

	// construct load vector
	double* F = new double[P];
	
	// for(int p=0; p<P; p++){
	// 	for(int q=0; q<Q; q++){
	// 		cout<<U[p][q]<<" ";
	// 	}
	// 	cout<<endl;
	// }

	write_matrix(P,Q,U);
	
	// cleanup
	for(int i = 0; i < Q; ++i) {
		delete [] U[i];
	}
	delete [] U;

	return 0;
}
