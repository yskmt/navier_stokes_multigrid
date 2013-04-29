#include "IO.h"

// write out the sparse matrix
int write_matrix(cuint P,
				 cuint Q,
				 double** U,
				 char* file_name)
{
	ofstream file_out;
	file_out.open (file_name);

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
				  double* F,
				  char* file_name)
{
	ofstream file_out;
	file_out.open (file_name);

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
				   cdouble dz,
				   cuint level
				   )
{
	char file_name[100];
	// write out the results in vtk format
	ofstream file_out;
	sprintf(file_name, "results_%i.vtk", level);
	file_out.open (file_name);
	if(!file_out.is_open()){
		return 1;
	}
	
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

	// output files for matlab
	// point coordinates and scalar result
	sprintf(file_name, "results_%i.dat", level);
	file_out.open(file_name);
	for(int n=0; n<n_dof; n++){
		one_d_to_three_d( n, I, J, i, j, k);
		file_out<<i*dx<<" "<<j*dy<<" "<<k*dz<<" "<<u[n]<<endl;
	}
	file_out.close();
	
}
