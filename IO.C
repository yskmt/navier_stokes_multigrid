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
int write_results(  double* U,
					double* V,
					double* W,
					double* P,
					cuint n_dof,
					cuint nx,
					cuint ny,
					cuint nz,
					cdouble xmin,
					cdouble ymin,
					cdouble zmin,
					cdouble hx,
					cdouble hy,
					cdouble hz,
					cuint ts,
					cdouble bcs[][6]
					)
{
	
	// build matrix with boundary conditions
	double* Ue = new double[(nx+1)*(ny+2)*(nz+2)];
	double* Ve = new double[(nx+2)*(ny+1)*(nz+2)];
	double* We = new double[(nx+2)*(ny+2)*(nz+1)];
	
	grid_matrix(U, Ue, nx-1, ny, nz, nx+1, ny+2, nz+2, X_DIR, bcs[0]);
	grid_matrix(V, Ve, nx, ny-1, nz, nx+2, ny+1, nz+2, Y_DIR, bcs[1]);
	grid_matrix(W, We, nx, ny, nz-1, nx+2, ny+2, nz+1, Z_DIR, bcs[2]);

	// get U, V, W defined at cell centers
	double* Uc = new double[(nx)*(ny+2)*(nz+2)];
	double* Vc = new double[(nx+2)*(ny)*(nz+2)];
	double* Wc = new double[(nx+2)*(ny+2)*(nz)];
	// average into cell centers
	average(Ue, Uc, nx+1, ny+2, nz+2, nx, ny+2, nz+2, X_DIR);
	average(Ve, Vc, nx+2, ny+1, nz+2, nx+2, ny, nz+2, Y_DIR);
	average(We, Wc, nx+2, ny+2, nz+1, nx+2, ny+2, nz, Z_DIR);
	
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
	
	// unsigned int i,j,k;
	// for(int n=0; n<n_dof; n++){
	// 	one_d_to_three_d( n, nx, ny, i, j, k);
	// 	file_out<<i*hx<<" "<<j*hy<<" "<<k*hz<<endl;
	// }

	// base for the grid
	cdouble xbase=xmin+hx/2;
	cdouble ybase=ymin+hy/2;
	cdouble zbase=zmin+hz/2;
	
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			for(int k=0; k<nz; k++){
				file_out<<xbase+i*hx<<" "<<ybase+j*hy<<" "<<zbase+k*hz<<endl;
			}
		}
	}
	
	file_out<<"POINT_DATA "<<n_dof<<endl;
	file_out<<"SCALARS P float 1"<<endl
			<<"LOOKUP_TABLE default"<<endl;
	// for(int n=0; n<n_dof; n++){
	// 	file_out<<P[n]<<endl;
	// }


	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			for(int k=0; k<nz; k++){
				uint t4;
				three_d_to_one_d(i,j,k, nx,ny, t4);
	
				file_out<<P[t4]<<endl;
			}
		}
	}

file_out<<"VECTORS velocity float"<<endl;

 for(int i=0; i<nx; i++){
	 for(int j=0; j<ny; j++){
		 for(int k=0; k<nz; k++){
			 uint t1, t2, t3;
			 three_d_to_one_d(i,j+1,k+1, nx, ny+2,  t1);
			 three_d_to_one_d(i+1,j,k+1, nx+2, ny,  t2);
			 three_d_to_one_d(i+1,j+1,k, nx+2, ny+2,t3);

			 file_out<<Uc[t1]<<" "
					 <<Vc[t2]<<" "
					 <<Wc[t3]<<endl;
		 }
	 }
 }


 // for(int n=0; n<n_dof; n++){
 // 		// uint i,j,k;
 // 		// one_d_to_three_d(n, nx, ny, i,j,k);
 // 		file_out<<Uc[n]<<" "
 // 				<<Vc[n]<<" "
 // 				<<Wc[n]<<endl;
 // 	}
	
	file_out.close();

	// output files for matlab
	// point coordinates and scalar result
	sprintf(file_name, "results_%i.dat", ts);
	file_out.open(file_name);

	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			for(int k=0; k<nz; k++){
				uint t1, t2, t3, t4;
				three_d_to_one_d(i,j+1,k+1, nx, ny+2,  t1);
				three_d_to_one_d(i+1,j,k+1, nx+2, ny,  t2);
				three_d_to_one_d(i+1,j+1,k, nx+2, ny+2,t3);
				three_d_to_one_d(i,j,k, nx,ny, t4);
				
				file_out<<xbase+i*hx<<" "<<ybase+j*hy<<" "<<zbase+k*hz<<" "
						<<P[t4]<<" "
						<<Uc[t1]<<" "
						<<Vc[t2]<<" "
						<<Wc[t3]<<endl;
			}
		}
	}
	
	// for(int n=0; n<n_dof; n++){
	// 	one_d_to_three_d( n, nx, ny, i, j, k);
	// 	file_out<<i*hx<<" "<<j*hy<<" "<<k*hz<<" "
	// 			<<P[n]<<" "
	// 			<<Uc[n]<<" "
	// 			<<Vc[n]<<" "
	// 			<<Wc[n]<<endl;
	// }
	file_out.close();
	
}

// write out 3d data for debuggin purpose
int write_3d_data( double* U,
				   cuint nx, cuint ny, cuint nz,
				   char* file_name )
{
	// output files for matlab
	ofstream file_out;
	file_out.open(file_name);

	if(!file_out.is_open()){
		cout<<"failed to open file"<<endl;
		return 1;
	}
	
	// 1d index
	uint t;

	// set up U
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			for(int k=0; k<nz; k++){
				three_d_to_one_d(i,j,k, nx, ny, t);
				file_out<<i<<" "<<j<<" "<<k<<" "<<U[t]<<endl;
			}
		}
	}

	file_out.close();
		
	return 0;
}
