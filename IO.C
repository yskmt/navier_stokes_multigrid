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


// write out the results
int write_results(  boost::multi_array<double, 3>& U,
				 boost::multi_array<double, 3>& V,
				 boost::multi_array<double, 3>& W,
				 double* P,
					cuint n_dof,
					cuint nx,
					cuint ny,
					cuint nz,
					cdouble hx,
					cdouble hy,
					cdouble hz,
					cuint ts,
					cdouble bcs[][6]
					)
{
	
	// build matrix with boundary conditions
	boost::multi_array<double, 3> Ue(boost::extents[nx+1][ny+2][nz+2]);
	boost::multi_array<double, 3> Ve(boost::extents[nx+2][ny+1][nz+2]);
	boost::multi_array<double, 3> We(boost::extents[nx+2][ny+2][nz+1]);
	grid_matrix(&U, &Ue, nx, ny, nz, X_DIR,
				bcs[0]);
	// grid_matrix(&V, &Ve, nx, ny, nz, Y_DIR,	bcs[1]);
	// grid_matrix(&W, &We, nx, ny, nz, Z_DIR, bcs[2]);

	// get U, V, W defined at cell centers
	boost::multi_array<double, 3> Uc(boost::extents[nx][ny+2][nz+2]);
	boost::multi_array<double, 3> Vc(boost::extents[nx+2][ny][nz+2]);
	boost::multi_array<double, 3> Wc(boost::extents[nx+2][ny+2][nz]);
	// average into cell centers
	average(Ue, Uc, X_DIR);
	average(Ve, Vc, Y_DIR);
	average(We, Wc, Z_DIR);
	

	char file_name[100];
	// write out the results in vtk format
	ofstream file_out;
	sprintf(file_name, "results_%i.vtk", ts);
	file_out.open (file_name);
	if(!file_out.is_open()){
		return 1;
	}
	
	// header
	file_out<<"# vtk DataFile Version 3.0"<<endl;
	file_out<<"3d incompressible NS problem"<<endl
			<<"ASCII"<<endl
			<<"DATASET STRUCTURED_GRID"<<endl
			<<"DIMENSIONS "<<nx<<" "<<ny<<" "<<nz<<endl
			<<"POINTS "<<n_dof<<" "<<"float"<<endl;
	
	unsigned int i,j,k;
	for(int n=0; n<n_dof; n++){
		one_d_to_three_d( n, nx, ny, i, j, k);
		file_out<<i*hx<<" "<<j*hy<<" "<<k*hz<<endl;
	}
	
	file_out<<"POINT_DATA "<<n_dof<<endl;
	file_out<<"SCALARS P float 1"<<endl
			<<"LOOKUP_TABLE default"<<endl;
	for(int n=0; n<n_dof; n++){
		file_out<<P[n]<<endl;
	}

	file_out<<"VECTORS U float"<<endl;
	for(int n=0; n<n_dof; n++){
		uint i,j,k;
		one_d_to_three_d(n, nx, ny, i,j,k);
		file_out<<Uc[i][j+1][k+1]<<" "
				<<Vc[i+1][j][k+1]<<" "
				<<Wc[i+1][j+1][k]<<endl;
	}
	
	file_out.close();

	// output files for matlab
	// point coordinates and scalar result
	sprintf(file_name, "results_%i.dat", ts);
	file_out.open(file_name);
	for(int n=0; n<n_dof; n++){
		one_d_to_three_d( n, nx, ny, i, j, k);
		file_out<<i*hx<<" "<<j*hy<<" "<<k*hz<<" "<<P[n]
				<<Uc[i][j+1][k+1]<<" "
				<<Vc[i+1][j][k+1]<<" "
				<<Wc[i+1][j+1][k]<<endl;
	}
	file_out.close();
	
}

// write out 3d data for debuggin purpose
int write_3d_data( boost::multi_array<double, 3>& U,
				  char* file_name )
{
	boost::multi_array_types::size_type const* sizes = U.shape();
	cuint nx = sizes[0];
	cuint ny = sizes[1];
	cuint nz = sizes[2];

	// output files for matlab
	ofstream file_out;
	file_out.open(file_name);

	// set up U
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			for(int k=0; k<nz; k++){
				file_out<<i<<" "<<j<<" "<<k<<" "<<U[i][j][k]<<endl;
			}
		}
	}

	file_out.close();
		
	return 0;
}
