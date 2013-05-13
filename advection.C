#include "advection.h"
#include "IO.h"

// treat nonlinear advection terms
void advection( double* U,
				double* V,
				double* W,
				cuint nx, cuint ny, cuint nz,
				cdouble hx, cdouble hy, cdouble hz,
				cdouble dt,
				cdouble bcs[][6]
				)
{	
	// upwinding parameter
	// cdouble gamma =
	// 	min( 1.2*dt*max(max(max_3d_array(U), max_3d_array(V)), max_3d_array(W)),
	// 		 1.0);

	// build matrix with boundary conditions
	double* Ue = new double[(nx+1)*(ny+2)*(nz+2)];
	double* Ve = new double[(nx+2)*(ny+1)*(nz+2)];
	double* We = new double[(nx+2)*(ny+2)*(nz+1)];
	grid_matrix(U, Ue, nx-1, ny, nz, nx+1, ny+2, nz+2, X_DIR, bcs[0]);
	grid_matrix(V, Ve, nx, ny-1, nz, nx+2, ny+1, nz+2, Y_DIR, bcs[1]);
	grid_matrix(W, We, nx, ny, nz-1, nx+2, ny+2, nz+1, Z_DIR, bcs[2]);

	// char Ufile [] = "U.dat";
	// write_3d_data(Ue, nx+1, ny+2, nz+2, Ufile);
	// char Vfile [] = "V.dat";
	// write_3d_data(Ve, nx+2, ny+1, nz+2, Vfile);
	// char Wfile [] = "W.dat";
	// write_3d_data(We, nx+2, ny+2, nz+1, Wfile);
	
	// calculate UV, UW, VW defined at mid-edges
	double* UV = new double[(nx+1)*(ny+1)*(nz+2)];
	double* UW = new double[(nx+1)*(ny+2)*(nz+1)];
	double* VW = new double[(nx+2)*(ny+1)*(nz+1)];
	calculate_edge_values(Ue, Ve, We, UV, UW, VW, nx, ny, nz);
		
	// calculate UV_y, UW_z, VU_x, VW_z, WU_x, WV_y
	// defined at corresponding points
	double*  UV_y = new double[(nx+1)*(ny)*(nz+2)];
	double*  UW_z = new double[(nx+1)*(ny+2)*(nz)];
	double*  VU_x = new double[(nx)*(ny+1)*(nz+2)];
	double*  VW_z = new double[(nx+2)*(ny+1)*(nz)];
	double*  WU_x = new double[(nx)*(ny+2)*(nz+1)];
	double*  WV_y = new double[(nx+2)*(ny)*(nz+1)];
	staggered_first_difference( UV, UV_y, nx+1, ny+1, nz+2,
								nx+1, ny, nz+2, hy, Y_DIR );
	staggered_first_difference( UW, UW_z, nx+1, ny+2, nz+1,
								nx+1, ny+2, nz, hz, Z_DIR );
	staggered_first_difference( UV, VU_x, nx+1, ny+1, nz+2,
								nx, ny+1, nz+2, hx, X_DIR );
	staggered_first_difference( VW, VW_z, nx+2, ny+1, nz+1,
								nx+2, ny+1, nz, hz, Z_DIR );
	staggered_first_difference( UW, WU_x, nx+1, ny+2, nz+1,
								nx, ny+2, nz+1, hx, X_DIR );
	staggered_first_difference( VW, WV_y, nx+2, ny+1, nz+1,
								nx+2, ny, nz+1, hy, Y_DIR );

	// get U, V, W defined at cell centers
	double* U2c = new double[(nx)*(ny+2)*(nz+2)];
	double* V2c = new double[(nx+2)*(ny)*(nz+2)];
	double* W2c = new double[(nx+2)*(ny+2)*(nz)];
	// average into cell centers
	average(Ue, U2c, nx+1, ny+2, nz+2, nx, ny+2, nz+2, X2_DIR);
	average(Ve, V2c, nx+2, ny+1, nz+2, nx+2, ny, nz+2, Y2_DIR);
	average(We, W2c, nx+2, ny+2, nz+1, nx+2, ny+2, nz, Z2_DIR);

	// get U2, V2, W2 defined at corresponding points
	double* U2_x = new double[(nx-1)*(ny+2)*(nz+2)];
	double* V2_y = new double[(nx+2)*(ny-1)*(nz+2)];
	double* W2_z = new double[(nx+2)*(ny+2)*(nz-1)];
	staggered_first_difference( U2c, U2_x, nx, ny+2, nz+2,
								nx-1, ny+2, nz+2, hx, X_DIR );
	staggered_first_difference( V2c, V2_y, nx+2, ny, nz+2,
								nx+2, ny-1, nz+2, hy, Y_DIR );
	staggered_first_difference( W2c, W2_z, nx+2, ny+2, nz,
								nx+2, ny+2, nz-1, hz, Z_DIR );
	
	// consolidate advection terms
	consolidate_advection( U, V, W,
						   U2_x, V2_y, W2_z,
						   UV_y, UW_z, VU_x, VW_z, WU_x, WV_y,
						   nx, ny, nz,
						   dt );

	// cleanup
	delete[] Ue, Ve, We;
	delete[] UV, UW, VW;
	delete[] UV_y, UW_z, VU_x, VW_z, WU_x, WV_y;
	delete[] U2c, V2c, W2c;
	delete[] U2_x, V2_y, W2_z;
	
	return;
}

// generate grid matrix
// dir: direction of velocity: 0:x 1:y 2:z
void grid_matrix( double* U,
				  double* Ue,
				  cuint nx, cuint ny, cuint nz,
				  cuint nxe, cuint nye, cuint nze,
				  cuint dir,
				  cdouble* bc )
{
	// 1d index
	uint t, t0, t1;

	// insert values
#pragma omp parallel for private(t, t0, t1) shared(Ue, U) num_threads(nt)
	for(int i=0; i<nxe; i++){
		for(int j=0; j<nye; j++){
			for(int k=0; k<nze; k++){
				three_d_to_one_d(i,j,k, nxe,nye, t);					
				if(i==0 || i==nxe-1 || j==0 || j==nye-1
				   || k==0 || k==nze-1)
					Ue[t] = 0.0;
				else{
					three_d_to_one_d(i-1,j-1,k-1, nx, ny, t1);
					Ue[t] = U[t1];
				}
			}
		}
	}

	// cout<<"U"<<endl;
	// if(dir==Y_DIR)
	// for(int i=0; i<(nx)*(ny)*(nz); i++)
	// 	cout<<U[i]<<endl;
		
	// for cofficients u
	if( dir==X_DIR){
		
#pragma omp parallel for private(t) shared(Ue) num_threads(nt)
		for(int j=0; j<nye; j++){
			for(int k=0; k<nze; k++){
				three_d_to_one_d(0,j,k, nxe,nye, t);
				Ue[t] = bc[0];
				three_d_to_one_d(nxe-1,j,k, nxe,nye, t);
				Ue[t] = bc[1];
			}
		}
		
		// y0, yl
#pragma omp parallel for private(t, t0, t1) shared(Ue) num_threads(nt)
		for(int i=1; i<nxe-1; i++){
			for(int k=1; k<nze-1; k++){
				three_d_to_one_d(i,0,k, nxe,nye, t0);
				three_d_to_one_d(i,1,k, nxe,nye, t1);
				Ue[t0] = 2*bc[2]-Ue[t1];
				
				three_d_to_one_d(i,nye-1,k, nxe,nye, t0);
				three_d_to_one_d(i,nye-2,k, nxe,nye, t1);
				Ue[t0] = 2*bc[3]-Ue[t1];
			}
		}

		// z0, zl
#pragma omp parallel for private(t, t0, t1) shared(Ue) num_threads(nt)
		for(int i=0; i<nxe; i++){
			for(int j=0; j<nye; j++){
				three_d_to_one_d(i,j,0, nxe,nye, t0);
				three_d_to_one_d(i,j,1, nxe,nye, t1);
				Ue[t0] = 2*bc[4]-Ue[t1];

				three_d_to_one_d(i,j,nze-1, nxe,nye, t0);
				three_d_to_one_d(i,j,nze-2, nxe,nye, t1);
				Ue[t0] = 2*bc[5]-Ue[t1];
			}
		}
		
	}

	// for cofficients v
	else if( dir==Y_DIR ){
		
		// account for boundary conditions
		// x0, xl
#pragma omp parallel for private(t, t0, t1) shared(Ue) num_threads(nt)
		for(int j=0; j<nye; j++){
			for(int k=0; k<nze; k++){
				three_d_to_one_d(0,j,k, nxe,nye, t0);
				three_d_to_one_d(1,j,k, nxe,nye, t1);
				Ue[t0] = 2*bc[0]-Ue[t1];
				
				three_d_to_one_d(nxe-1,j,k, nxe,nye, t0);
				three_d_to_one_d(nxe-2,j,k, nxe,nye, t1);
				Ue[t0] = 2*bc[0]-Ue[t1];

			}
		}

		// y0, yl
#pragma omp parallel for private(t) shared(Ue) num_threads(nt)
		for(int i=0; i<nxe; i++){
			for(int k=0; k<nze; k++){
				three_d_to_one_d(i,0,k, nxe,nye, t);
				Ue[t] = bc[2];
				three_d_to_one_d(i,nye-1,k, nxe,nye, t);
				Ue[t] = bc[3];

			}
		}

		// z0, zl
#pragma omp parallel for private(t, t0, t1) shared(Ue) num_threads(nt)
		for(int i=0; i<nxe; i++){
			for(int j=0; j<nye; j++){
				three_d_to_one_d(i,j,0, nxe,nye, t0);
				three_d_to_one_d(i,j,1, nxe,nye, t1);
				Ue[t0] = 2*bc[4]-Ue[t1];
		
				three_d_to_one_d(i,j,nze-1, nxe,nye, t0);
				three_d_to_one_d(i,j,nze-2, nxe,nye, t1);
				Ue[t0] = 2*bc[5]-Ue[t1];
		
			}
		}
		
	}

	// for cofficients w
	else if( dir==Z_DIR ){

		// account for boundary conditions
		// x0, xl
#pragma omp parallel for private(t, t0, t1) shared(Ue) num_threads(nt)
		for(int j=0; j<nye; j++){
			for(int k=0; k<nze; k++){
				three_d_to_one_d(0,j,k, nxe,nye, t0);
				three_d_to_one_d(1,j,k, nxe,nye, t1);
				Ue[t0] = 2*bc[0]-Ue[t1];

				three_d_to_one_d(nxe-1,j,k, nxe,nye, t0);
				three_d_to_one_d(nxe-2,j,k, nxe,nye, t1);
				Ue[t0] = 2*bc[1]-Ue[t1];
			}
		}

		// y0, yl
#pragma omp parallel for private(t, t0, t1) shared(Ue) num_threads(nt)
		for(int i=0; i<nxe; i++){
			for(int k=0; k<nze; k++){
				three_d_to_one_d(i,0,k, nxe,nye, t0);
				three_d_to_one_d(i,1,k, nxe,nye, t1);
				Ue[t0] = 2*bc[2]-Ue[t1];
				
				three_d_to_one_d(i,nye-1,k, nxe,nye, t0);
				three_d_to_one_d(i,nye-2,k, nxe,nye, t1);
				Ue[t0] = 2*bc[3]-Ue[t1];
			}
		}

		// z0, zl
#pragma omp parallel for private(t) shared(Ue) num_threads(nt)
		for(int i=0; i<nxe; i++){
			for(int j=0; j<nye; j++){
				three_d_to_one_d(i,j,0, nxe,nye, t);
				Ue[t] = bc[4];
				three_d_to_one_d(i,j,nze-1, nxe,nye, t);
				Ue[t] = bc[5];
			}
		}
	}
	
	return;	
}

// average values in whatever way
// dir: averaging direction:
// 0: x-direction
// 1: y-direction
// 2: z-direction
// 3: xy-direction
// 4: xz-direction
// 5: yz-direction
// 6: x-direction and square
// 7: y-direction and square
// 8: z-direction and square
void average( const double* Ue, // raw value
			  double* Ua, // averaged value
			  cuint nxe, cuint nye, cuint nze,
			  cuint nxa, cuint nya, cuint nza,
			  cuint dir )
{
	uint ta, te1, te2;

	// x-average
	if(dir==X_DIR){
		// interpolate values by averaging
#pragma omp parallel for private(ta, te1, te2) shared(Ue, Ua) num_threads(nt)
		for(int i=0; i<(nxa); i++){
			for(int j=0; j<(nya); j++){
				for(int k=0; k<(nza); k++){
					three_d_to_one_d(i,j,k, nxa, nya, ta);
					three_d_to_one_d(i+1,j,k, nxe, nye, te1);
					three_d_to_one_d(i,j,k, nxe, nye, te2);
					Ua[ta] = (Ue[te1]+Ue[te2])/2;
				}
			}
		}
	}

	// y-average
	else if(dir==Y_DIR){
		// interpolate values by averaging
#pragma omp parallel for private(ta, te1, te2) shared(Ue, Ua) num_threads(nt)
		for(int i=0; i<(nxa); i++){
			for(int j=0; j<(nya); j++){
				for(int k=0; k<(nza); k++){
					three_d_to_one_d(i,j,k, nxa, nya, ta);
					three_d_to_one_d(i,j+1,k, nxe, nye, te1);
					three_d_to_one_d(i,j,k, nxe, nye, te2);
					Ua[ta] = (Ue[te1]+Ue[te2])/2;
				}
			}
		}
	}

	// z-average
	else if(dir==Z_DIR){
		// interpolate values by averaging
#pragma omp parallel for private(ta, te1, te2) shared(Ue, Ua) num_threads(nt)
		for(int i=0; i<(nxa); i++){
			for(int j=0; j<(nya); j++){
				for(int k=0; k<(nza); k++){
					three_d_to_one_d(i,j,k, nxa, nya, ta);
					three_d_to_one_d(i,j,k+1, nxe, nye, te1);
					three_d_to_one_d(i,j,k, nxe, nye, te2);
					Ua[ta] = (Ue[te1]+Ue[te2])/2;
				}
			}
		}
	}

	// // xy-direction
	// else if(dir==3){
	// 	// interpolate values by averaging
	// 	for(int i=0; i<(nx-1); i++){
	// 		for(int j=0; j<(ny-1); j++){
	// 			for(int k=0; k<(nz); k++){
	// 				Ua[i][j][k] = (Ue[i][j][k]+Ue[i+1][j][k]
	// 							   + Ue[i][j+1][k]+Ue[i+1][j+1][k])/4;
	// 			}
	// 		}
	// 	}
	// }
	
	// // xz-direction
	// else if(dir==4){
	// 	// interpolate values by averaging
	// 	for(int i=0; i<(nx-1); i++){
	// 		for(int j=0; j<(ny); j++){
	// 			for(int k=0; k<(nz-1); k++){
	// 				Ua[i][j][k] = (Ue[i][j][k]+Ue[i+1][j][k]
	// 							   + Ue[i][j][k+1]+Ue[i+1][j][k+1])/4;
	// 			}
	// 		}
	// 	}
	// }

	// // yz-direction
	// else if(dir==5){
	// 	// interpolate values by averaging
	// 	for(int i=0; i<(nx); i++){
	// 		for(int j=0; j<(ny-1); j++){
	// 			for(int k=0; k<(nz-1); k++){
	// 				Ua[i][j][k] = (Ue[i][j][k]+Ue[i][j+1][k]
	// 							   + Ue[i][j][k+1]+Ue[i][j+1][k+1])/4;
	// 			}
	// 		}
	// 	}
	// }

	// x-direction and square
	else if(dir==X2_DIR){
		// interpolate values by averaging
#pragma omp parallel for private(ta, te1, te2) shared(Ue, Ua) num_threads(nt)
		for(int i=0; i<(nxa); i++){
			for(int j=0; j<(nya); j++){
				for(int k=0; k<(nza); k++){
					three_d_to_one_d(i,j,k, nxa, nya, ta);
					three_d_to_one_d(i+1,j,k, nxe, nye, te1);
					three_d_to_one_d(i,j,k, nxe, nye, te2);
					
					Ua[ta] = pow((Ue[te1]+Ue[te2])/2,2);
				}
			}
		}
	}

	// y-direction and square
	else if(dir==Y2_DIR){
		// interpolate values by averaging
#pragma omp parallel for private(ta, te1, te2) shared(Ue, Ua) num_threads(nt)
		for(int i=0; i<(nxa); i++){
			for(int j=0; j<(nya); j++){
				for(int k=0; k<(nza); k++){
					three_d_to_one_d(i,j,k, nxa, nya, ta);
					three_d_to_one_d(i,j+1,k, nxe, nye, te1);
					three_d_to_one_d(i,j,k, nxe, nye, te2);
					
					Ua[ta] = pow((Ue[te1]+Ue[te2])/2,2);
				}
			}
		}
	}

	// z-direction and square
	if(dir==Z2_DIR){
		// interpolate values by averaging
#pragma omp parallel for private(ta, te1, te2) shared(Ue, Ua) num_threads(nt)
		for(int i=0; i<(nxa); i++){
			for(int j=0; j<(nya); j++){
				for(int k=0; k<(nza); k++){
					three_d_to_one_d(i,j,k, nxa, nya, ta);
					three_d_to_one_d(i,j,k+1, nxe, nye, te1);
					three_d_to_one_d(i,j,k, nxe, nye, te2);
					
					Ua[ta] = pow((Ue[te1]+Ue[te2])/2,2);
				}
			}
		}
	}

	
	return;
}

// get maximum of the abs values of 3d array
double max_3d_array( const boost::multi_array<double, 3>& U )
{
	boost::multi_array_types::size_type const* sizes = U.shape();
	cuint nx = sizes[0];
	cuint ny = sizes[1];
	cuint nz = sizes[2];
	
	double max_value = abs(U[0][0][0]);

	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			for(int k=0; k<nz; k++){
				if(max_value < abs(U[i][j][k]))
					max_value = abs(U[i][j][k]);
			}
		}
	}

	return max_value;
}

// get 1d staggered difference
// Ua is an averaged U value at the cell vertices
// get the difference value at the center of cell
void staggered_first_difference( const double* UV,
								 double* UV_x,
								 cuint nx, cuint ny, cuint nz,
								 cuint nx_x, cuint ny_x, cuint nz_x,
								 cdouble h,
								 cuint dir
								 )
{
	uint t1, t2, t3;
	// difference in x-direction
	if( dir==X_DIR){
#pragma omp parallel for private(t1, t2, t3) shared(UV_x, UV) num_threads(nt)
		for(int i=0; i<nx-1; i++){
			for(int j=0; j<ny; j++){
				for(int k=0; k<nz; k++){
					three_d_to_one_d(i,j,k, nx_x, ny_x, t1);
					three_d_to_one_d(i+1,j,k, nx, ny, t2);
					three_d_to_one_d(i,j,k, nx, ny, t3);
					
					UV_x[t1] = (UV[t2]-UV[t3])/h;
				}
			}
		}
	}

	// difference in y-direction
	if( dir==Y_DIR ){
#pragma omp parallel for private(t1, t2, t3) shared(UV_x, UV) num_threads(nt)
		for(int i=0; i<nx; i++){
			for(int j=0; j<ny-1; j++){
				for(int k=0; k<nz; k++){
					three_d_to_one_d(i,j,k, nx_x, ny_x, t1);
					three_d_to_one_d(i,j+1,k, nx, ny, t2);
					three_d_to_one_d(i,j,k, nx, ny, t3);
					
					UV_x[t1] = (UV[t2]-UV[t3])/h;
				}
			}
		}
	}

	// difference in z-direction
	if( dir==Z_DIR ){
#pragma omp parallel for private(t1, t2, t3) shared(UV_x, UV) num_threads(nt)
		for(int i=0; i<nx; i++){
			for(int j=0; j<ny; j++){
				for(int k=0; k<nz-1; k++){
					three_d_to_one_d(i,j,k, nx_x, ny_x, t1);
					three_d_to_one_d(i,j,k+1, nx, ny, t2);
					three_d_to_one_d(i,j,k, nx, ny, t3);
			
					UV_x[t1] = (UV[t2]-UV[t3])/h;
				}
			}
		}
	}
}
	

// get upwinding differences
void upwind_difference( const boost::multi_array<double, 3>& U,
						boost::multi_array<double, 3>& Ud,
						cuint dir )
{
	boost::multi_array_types::size_type const* sizes = U.shape();
	cuint nx = sizes[0];
	cuint ny = sizes[1];
	cuint nz = sizes[2];

	// differece in x-direction
	if(dir==X_DIR){
		for(int i=0; i<nx-1; i++){
			for(int j=0; j<ny; j++){
				for(int k=0; k<nz; k++){
					Ud[i][j][k] = (U[i+1][j][k]-U[i][j][k]) / 2;
				}
			}
		}
	}

   	// differece in y-direction
	else if(dir==Y_DIR){
		for(int i=0; i<nx; i++){
			for(int j=0; j<ny-1; j++){
				for(int k=0; k<nz; k++){
					Ud[i][j][k] = (U[i][j+1][k]-U[i][j][k]) / 2;
				}
			}
		}
	}
	// differece in z-direction
	else if (dir==Z_DIR){
		for(int i=0; i<nx; i++){
			for(int j=0; j<ny; j++){
				for(int k=0; k<nz-1; k++){
					Ud[i][j][k] = (U[i][j][k+1]-U[i][j][k]) / 2;
				}
			}
		}

	}
	
	return;
}

// get central first difference at center of element
void central_first_difference( const boost::multi_array<double, 3>& U2,
							   boost::multi_array<double, 3>& U2_x,
							   cdouble h,
							   cuint dir )
{
	boost::multi_array_types::size_type const* sizes = U2.shape();
	cuint nx = sizes[0];
	cuint ny = sizes[1];
	cuint nz = sizes[2];

	// x-difference
	if( dir==X_DIR){
		for(int i=0; i<nx-2; i++){
			for(int j=0; j<ny; j++){
				for(int k=0; k<nz; k++){	
					U2_x[i][j][k] = (U2[i][j][k]+U2[i+2][j][k])/(2*h);
				}
			}
		}
	}

	// y-difference
	if( dir==Y_DIR){
		for(int i=0; i<nx; i++){
			for(int j=0; j<ny-2; j++){
				for(int k=0; k<nz; k++){	
					U2_x[i][j][k] = (U2[i][j][k]+U2[i][j+2][k])/(2*h);
				}
			}
		}
	}

	// z-difference
	if( dir==Z_DIR){
		for(int i=0; i<nx; i++){
			for(int j=0; j<ny; j++){
				for(int k=0; k<nz-2; k++){	
					U2_x[i][j][k] = (U2[i][j][k]+U2[i][j][k+2])/(2*h);
				}
			}
		}
	}

	return;
}

// get mixed edge values
void calculate_edge_values( double* Ue,
							double* Ve,
							double* We,
							double* UV,
							double* UW,
							double* VW,
							cuint nx, cuint ny, cuint nz)
{
	// average each value
	double* Uay = new double[(nx+1)*(ny+1)*(nz+2)];
	double* Vax = new double[(nx+1)*(ny+1)*(nz+2)];
	double* Uaz = new double[(nx+1)*(ny+2)*(nz+1)];
	double* Wax = new double[(nx+1)*(ny+2)*(nz+1)];
	double* Vaz = new double[(nx+2)*(ny+1)*(nz+1)];
	double* Way = new double[(nx+2)*(ny+1)*(nz+1)];
	
	average(Ue, Uay, nx+1, ny+2, nz+2, nx+1, ny+1, nz+2, Y_DIR);
	average(Ve, Vax, nx+2, ny+1, nz+2, nx+1, ny+1, nz+2, X_DIR);
	average(Ue, Uaz, nx+1, ny+2, nz+2, nx+1, ny+2, nz+1, Z_DIR);
	average(We, Wax, nx+2, ny+2, nz+1, nx+1, ny+2, nz+1, X_DIR);
	average(Ve, Vaz, nx+2, ny+1, nz+2, nx+2, ny+1, nz+1, Z_DIR);
	average(We, Way, nx+2, ny+2, nz+1, nx+2, ny+1, nz+1, Y_DIR);

	uint t;
	
#pragma omp parallel for shared(UV, Uay, Vax) num_threads(nt)
	for(int i=0; i<nx+1; i++){
		for(int j=0; j<ny+1; j++){
			for(int k=0; k<nz+2; k++){
				three_d_to_one_d(i,j,k, nx+1, ny+1, t);
				UV[t] = Uay[t]*Vax[t];
			}
		}
	}

#pragma omp parallel for shared(UW, Uaz, Wax) num_threads(nt)
	for(int i=0; i<nx+1; i++){
		for(int j=0; j<ny+2; j++){
			for(int k=0; k<nz+1; k++){
				three_d_to_one_d(i,j,k, nx+1, ny+2, t);
				UW[t] = Uaz[t]*Wax[t];
			}
		}
	}

#pragma omp parallel for shared(VW, Vaz, Wax) num_threads(nt)
	for(int i=0; i<nx+2; i++){
		for(int j=0; j<ny+1; j++){
			for(int k=0; k<nz+1; k++){
				three_d_to_one_d(i,j,k, nx+2, ny+1, t);
				VW[t] = Vaz[t]*Way[t];
			}
		}
	}

	//cleanup
	delete[] Uay, Vax, Uaz, Wax, Vaz, Way;
	
	return;

}

// consolidate advection terms
void consolidate_advection( double* U,
							double* V,
							double* W,
							double* U2_x,
							double* V2_y,
							double* W2_z,
							double* UV_y,
							double* UW_z,
							double* VU_x,
							double* VW_z,
							double* WU_x,
							double* WV_y,
							cuint nx, cuint ny, cuint nz,
							cdouble dt )
{
	// need to truncate some terms
	uint t0, t1, t2, t3;

	
	// x-direction
#pragma omp parallel for private(t0, t1, t2, t3) shared(U, U2_x, UV_y, UW_z) num_threads(nt)
	for(int i=0; i<(nx-1); i++){
		for(int j=0; j<(ny); j++){
			for(int k=0; k<(nz); k++){
				three_d_to_one_d(i,j,k, nx-1,ny, t0);
				three_d_to_one_d(i,j+1,k+1, nx-1,ny+2, t1);
				three_d_to_one_d(i+1,j,k+1, nx+1,ny, t2);
				three_d_to_one_d(i+1,j+1,k, nx+1,ny+2, t3);
				
				U[t0] = U[t0] - dt * (U2_x[t1] + UV_y[t2] + UW_z[t3]);
			}
		}
	}

	// y-direction
#pragma omp parallel for private(t0, t1, t2, t3) shared(V, V2_y, VU_x, VW_z) num_threads(nt)
	for(int i=0; i<(nx); i++){
		for(int j=0; j<(ny-1); j++){
			for(int k=0; k<(nz); k++){
				three_d_to_one_d(i,j,k, nx,ny-1, t0);
				three_d_to_one_d(i+1,j,k+1, nx+2,ny-1, t1);
				three_d_to_one_d(i,j+1,k+1, nx,ny+1, t2);
				three_d_to_one_d(i+1,j+1,k, nx+2,ny+1, t3);
								
				V[t0] = V[t0] - dt * (V2_y[t1] + VU_x[t2] + VW_z[t3]);
			}
		}
	}

	// z-direction
#pragma omp parallel for private(t0, t1, t2, t3) shared(W, W2_z, WU_x, WV_y) num_threads(nt)
	for(int i=0; i<(nx); i++){
		for(int j=0; j<(ny); j++){
			for(int k=0; k<(nz-1); k++){
				three_d_to_one_d(i,j,k, nx,ny, t0);
				three_d_to_one_d(i+1,j+1,k, nx+2,ny+2, t1);
				three_d_to_one_d(i,j+1,k+1, nx, ny+2, t2);
				three_d_to_one_d(i+1,j,k+1, nx+2, ny, t3);
								
				W[t0] = W[t0] - dt *
					(W2_z[t1] + WU_x[t2] + WV_y[t3]);
			}
		}
	}
	
	return;

}
