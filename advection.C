#include "advection.h"

// set up initial conditions
void initial_conditoins()
{

	return;
	
}

// treat nonlinear advection terms
void advection( boost::multi_array<double, 3>& U,
				boost::multi_array<double, 3>& V,
				boost::multi_array<double, 3>& W,
				cuint nx, cuint ny, cuint nz,
				cdouble hx, cdouble hy, cdouble hz,
				cdouble dt
				)
{
	
	// boundary values
	double uB = 0.0; //bottom
	double uT = 0.0; //top

	// upwinding parameter
	cdouble gamma =
		min( 1.2*dt*max(max(max_3d_array(U), max_3d_array(V)), max_3d_array(W)),
			 1.0);

	// build matrix with boundary conditions
	boost::multi_array<double, 3> Ue(boost::extents[nx+1][ny+2][nz+2]);
	boost::multi_array<double, 3> Ve(boost::extents[nx+2][ny+1][nz+2]);
	boost::multi_array<double, 3> We(boost::extents[nx+2][ny+2][nz+1]);
	grid_matrix(Ue, nx, ny, nz, X_DIR, 0,0,0,0, uB, uT);
	grid_matrix(Ve, nx, ny, nz, Y_DIR, 0,0,0,0, uB, uT);
	grid_matrix(We, nx, ny, nz, Z_DIR, 0,0,0,0, uB, uT);
	
	// calculate UV, UW, VW defined at mid-edges
	boost::multi_array<double, 3> UV(boost::extents[nx+1][ny+1][nz+2]);
	boost::multi_array<double, 3> UW(boost::extents[nx+1][ny+2][nz+1]);
	boost::multi_array<double, 3> VW(boost::extents[nx+2][ny+1][nz+1]);
	calculate_edge_values(Ue, Ve, We, UV, UW, VW, nx, ny, nz);
		
	// calculate UV_y, UW_z, VU_x, VW_z, WU_x, WV_y
	// defined at corresponding points
	boost::multi_array<double, 3> UV_y(boost::extents[nx+1][ny][nz+2]);
	boost::multi_array<double, 3> UW_z(boost::extents[nx+1][ny+2][nz]);
	boost::multi_array<double, 3> VU_x(boost::extents[nx][ny+1][nz+2]);
	boost::multi_array<double, 3> VW_z(boost::extents[nx+2][ny+1][nz]);
	boost::multi_array<double, 3> WU_x(boost::extents[nx][ny+2][nz+1]);
	boost::multi_array<double, 3> WV_y(boost::extents[nx+2][ny][nz+1]);
	staggered_first_difference( UV, UV_y, hy, Y_DIR );
	staggered_first_difference( UW, UW_z, hz, Z_DIR );
	staggered_first_difference( UV, VU_x, hx, X_DIR );
	staggered_first_difference( VW, VW_z, hz, Z_DIR );
	staggered_first_difference( UW, WU_x, hx, X_DIR );
	staggered_first_difference( VW, WV_y, hy, Y_DIR );

	// get U, V, W defined at cell centers
	boost::multi_array<double, 3> U2c(boost::extents[nx][ny+2][nz+2]);
	boost::multi_array<double, 3> V2c(boost::extents[nx+2][ny][nz+2]);
	boost::multi_array<double, 3> W2c(boost::extents[nx+2][ny+2][nz]);
	// average into cell centers
	average(Ue, U2c, X2_DIR);
	average(Ve, V2c, Y2_DIR);
	average(We, W2c, Z2_DIR);

	// get U2, V2, W2 defined at corresponding points
	boost::multi_array<double, 3> U2_x(boost::extents[nx-1][ny+2][nz+2]);
	boost::multi_array<double, 3> V2_y(boost::extents[nx+2][ny-1][nz+2]);
	boost::multi_array<double, 3> W2_z(boost::extents[nx+2][ny+2][nz-1]);
	staggered_first_difference(  U2c, U2_x, hx, X_DIR );
	staggered_first_difference(  V2c, V2_y, hy, Y_DIR );
	staggered_first_difference(  W2c, W2_z, hz, Z_DIR );
	
	// consolidate advection terms
	consolidate_advection( U, V, W,
						   U2_x, V2_y, W2_z,
						   UV_y, UW_z, VU_x, VW_z, WU_x, WV_y,
						   nx, ny, nz,
						   dt );
	
	return;
}

// generate grid matrix
// dir: direction of velocity: 0:x 1:y 2:z
void grid_matrix( boost::multi_array<double, 3>& Ue,
				  cuint nx, cuint ny, cuint nz,
				  cuint dir,
				  cdouble u0, // x0
				  cdouble u1, // xl
				  cdouble u2, // y0
				  cdouble u3, // yl
				  cdouble u4, // z0
				  cdouble u5 ) //zl
	
{
	double x_size, y_size, z_size;
	
	// u velocity
	if(dir==0){
		x_size=nx+1; y_size=ny+2; z_size=nz+2;
	}
	// v velocity
	else if(dir==1){
		x_size=nx+2; y_size=ny+1; z_size=nz+2;
	}
	// w velocity
	else if(dir==2){
		x_size=nx+2; y_size=ny+2; z_size=nz+1;
	}
	else{
		cout<<"ERROR: direction "<< dir <<" does not exist."<<endl;
		return;
	}
	
	// for cofficients u
	if( dir==0){
		// account for boundary conditions
		// x0, xl
		for(int j=0; j<y_size; j++){
			for(int k; k<z_size; k++){
				Ue[0][j][k] = u0;
				Ue[x_size-1][j][k] = u1;
			}
		}
		// y0, yl
		for(int i=0; i<x_size; i++){
			for(int k; k<z_size; k++){
				Ue[i][0][k] = 2*u2-Ue[i][1][k];
				Ue[i][y_size-1][k] = 2*u3-Ue[i][y_size-2][k];
			
			}
		}
		// z0, zl
		for(int i=0; i<x_size; i++){
			for(int j; j<y_size; j++){
				Ue[i][j][0] = 2*u4-Ue[i][j][1];
				Ue[i][j][z_size-1] = 2*u5-Ue[i][j][z_size-2];
			}
		}
	}

	// for cofficients v
	else if( dir==1 ){
		// account for boundary conditions
		// x0, xl
		for(int j=0; j<y_size; j++){
			for(int k; k<z_size; k++){
				Ue[0][j][k] = 2*u0-Ue[1][j][k];
				Ue[x_size-1][j][k] = 2*u1-Ue[x_size-2][j][k];
			}
		}
		// y0, yl
		for(int i=0; i<x_size; i++){
			for(int k; k<z_size; k++){
				Ue[i][0][k] = u2;
				Ue[i][y_size-1][k] = u3;
			}
		}
		// z0, zl
		for(int i=0; i<x_size; i++){
			for(int j; j<y_size; j++){
				Ue[i][j][0] = 2*u4-Ue[i][j][1];
				Ue[i][j][z_size-1] = 2*u5-Ue[i][j][z_size-2];
			}
		}
	}

	// for cofficients w
	else if( dir==2 ){
		// account for boundary conditions
		// x0, xl
		for(int j=0; j<y_size; j++){
			for(int k; k<z_size; k++){
				Ue[0][j][k] = 2*u0-Ue[1][j][k];
				Ue[x_size-1][j][k] = 2*u1-Ue[x_size-2][j][k];
			}
		}
		// y0, yl
		for(int i=0; i<x_size; i++){
			for(int k; k<z_size; k++){
				Ue[i][0][k] = 2*u2-Ue[i][1][k];
				Ue[i][y_size-1][k] = 2*u3-Ue[i][y_size-2][k];
			}
		}
		// z0, zl
		for(int i=0; i<x_size; i++){
			for(int j; j<y_size; j++){
				Ue[i][j][0] = u4;
				Ue[i][j][z_size-1] = u5;
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
void average( const boost::multi_array<double, 3>& Ue,
			  boost::multi_array<double, 3>& Ua,
			  cuint dir )
{
	// boost::size_type Ua.shape();
	boost::multi_array_types::size_type const* sizes = Ua.shape();
	cuint nx = sizes[0];
	cuint ny = sizes[1];
	cuint nz = sizes[2];

	if(dir==0){
		// interpolate values by averaging
		for(int i=0; i<(nx-1); i++){
			for(int j=0; j<(ny); j++){
				for(int k=0; k<(nz); k++){
					Ua[i][j][k] = (Ue[i+1][j][k]+Ue[i][j][k])/2;
				}
			}
		}
	}

	else if(dir==1){
		// interpolate values by averaging
		for(int i=0; i<(nx); i++){
			for(int j=0; j<(ny-1); j++){
				for(int k=0; k<(nz); k++){
					Ua[i][j][k] = (Ue[i][j+1][k]+Ue[i][j][k])/2;
				}
			}
		}
	}
	
	else if(dir==2){
		// interpolate values by averaging
		for(int i=0; i<(nx); i++){
			for(int j=0; j<(ny); j++){
				for(int k=0; k<(nz-1); k++){
					Ua[i][j][k] = (Ue[i][j][k+1]+Ue[i][j][k])/2;
				}
			}
		}
	}

	// xy-direction
	else if(dir==3){
		// interpolate values by averaging
		for(int i=0; i<(nx-1); i++){
			for(int j=0; j<(ny-1); j++){
				for(int k=0; k<(nz); k++){
					Ua[i][j][k] = (Ue[i][j][k]+Ue[i+1][j][k]
								   + Ue[i][j+1][k]+Ue[i+1][j+1][k])/4;
				}
			}
		}
	}
	
	// xz-direction
	else if(dir==4){
		// interpolate values by averaging
		for(int i=0; i<(nx-1); i++){
			for(int j=0; j<(ny); j++){
				for(int k=0; k<(nz-1); k++){
					Ua[i][j][k] = (Ue[i][j][k]+Ue[i+1][j][k]
								   + Ue[i][j][k+1]+Ue[i+1][j][k+1])/4;
				}
			}
		}
	}

	// yz-direction
	else if(dir==5){
		// interpolate values by averaging
		for(int i=0; i<(nx); i++){
			for(int j=0; j<(ny-1); j++){
				for(int k=0; k<(nz-1); k++){
					Ua[i][j][k] = (Ue[i][j][k]+Ue[i][j+1][k]
								   + Ue[i][j][k+1]+Ue[i][j+1][k+1])/4;
				}
			}
		}
	}

	// x-direction and square
	else if(dir==6){
		// interpolate values by averaging
		for(int i=0; i<(nx-1); i++){
			for(int j=0; j<(ny); j++){
				for(int k=0; k<(nz); k++){
					Ua[i][j][k] = pow((Ue[i+1][j][k]+Ue[i][j][k])/2,2);
				}
			}
		}
	}

	// y-direction and square
	else if(dir==7){
		// interpolate values by averaging
		for(int i=0; i<(nx); i++){
			for(int j=0; j<(ny-1); j++){
				for(int k=0; k<(nz); k++){
					Ua[i][j][k] = pow((Ue[i][j+1][k]+Ue[i][j][k])/2,2);
				}
			}
		}
	}

	// z-direction and square
	if(dir==8){
		// interpolate values by averaging
		for(int i=0; i<(nx); i++){
			for(int j=0; j<(ny); j++){
				for(int k=0; k<(nz-1); k++){
					Ua[i][j][k] = pow((Ue[i][j][k+1]+Ue[i][j][k])/2,2);
				}
			}
		}
	}

	
	return;
}

// first order 
void diff_1()
{


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
void staggered_first_difference( const boost::multi_array<double, 3>& UV,
								 boost::multi_array<double, 3>& UV_x,
								 cdouble h,
								 cuint dir
								 )
{
	boost::multi_array_types::size_type const* sizes = UV.shape();
	cuint nx = sizes[0];
	cuint ny = sizes[1];
	cuint nz = sizes[2];
	
	// difference in x-direction
	if( dir==X_DIR){
		// boost::multi_array<double, 3> UV_x(boost::extents[nx-1][ny][nz]);
		for(int i=0; i<nx-1; i++){
			for(int j=0; j<ny; j++){
				for(int k=0; k<nz; k++){
					UV_x[i][j][k] = (UV[i+1][j][k]-UV[i][j][k])/h;
				}
			}
		}
	}

	// difference in y-direction
	if( dir==Y_DIR ){
		// boost::multi_array<double, 3> UV_y(boost::extents[nx][ny-1][nz]);
		for(int i=0; i<nx; i++){
			for(int j=0; j<ny-1; j++){
				for(int k=0; k<nz; k++){
					UV_x[i][j][k] = (UV[i][j+1][k]-UV[i][j][k])/h;
				}
			}
		}
	}

	// difference in z-direction
	if( dir==Z_DIR ){
		// boost::multi_array<double, 3> UV_z(boost::extents[nx][ny][nz-1]);
		for(int i=0; i<nx; i++){
			for(int j=0; j<ny; j++){
				for(int k=0; k<nz-1; k++){
					UV_x[i][j][k] = (UV[i][j][k+1]-UV[i][j][k])/h;
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
void calculate_edge_values( const boost::multi_array<double, 3>& Ue,
							const boost::multi_array<double, 3>& Ve,
							const boost::multi_array<double, 3>& We,
							boost::multi_array<double, 3>& UV,
							boost::multi_array<double, 3>& UW,
							boost::multi_array<double, 3>& VW,
							cuint nx, cuint ny, cuint nz)
{
	boost::multi_array<double, 3> Uay(boost::extents[nx+1][ny+1][nz+2]);
	boost::multi_array<double, 3> Vax(boost::extents[nx+1][ny+1][nz+2]);
	boost::multi_array<double, 3> Uaz(boost::extents[nx+1][ny+2][nz+1]);
	boost::multi_array<double, 3> Wax(boost::extents[nx+1][ny+2][nz+1]);
	average(Ue, Uay, Y_DIR);
	average(Ve, Vax, X_DIR);
	average(Ue, Uaz, Z_DIR);
	average(Ue, Wax, X_DIR);

	// boost::multi_array<double, 3> Va_x(boost::extents[nx+1][ny+1][nz+2]);
	// boost::multi_array<double, 3> Ua_y(boost::extents[nx+1][ny+1][nz+2]);
	boost::multi_array<double, 3> Vaz(boost::extents[nx+2][ny+1][nz+1]);
	boost::multi_array<double, 3> Way(boost::extents[nx+2][ny+1][nz+1]);
   	// average(Ve, Va_x, X_DIR);
	// average(Ue, Ua_y, Y_DIR);
	average(Ve, Vaz, Z_DIR);
	average(We, Way, Y_DIR);
	
	// boost::multi_array<double, 3> Wa_x(boost::extents[nx+1][ny+2][nz+1]);
	// boost::multi_array<double, 3> Ua_z(boost::extents[nx+1][ny+2][nz+1]);
	// boost::multi_array<double, 3> Wa_y(boost::extents[nx+2][ny+1][nz+1]);
	// boost::multi_array<double, 3> Va_z(boost::extents[nx+2][ny+1][nz+1]);
	// average(We, Wa_x, X_DIR);
	// average(Ue, Ua_z, Z_DIR);
	// average(We, Wa_y, Y_DIR);
	// average(Ve, Va_z, Z_DIR);
	
	for(int i=0; i<nx+1; i++){
		for(int j=0; j<ny+1; j++){
			for(int k=0; k<nz+2; k++){
				UV[i][j][k] = Uay[i][j][k]*Vax[i][j][k];
			}
		}
	}

	for(int i=0; i<nx+1; i++){
		for(int j=0; j<ny+2; j++){
			for(int k=0; k<nz+1; k++){
				UW[i][j][k] = Uaz[i][j][k]*Wax[i][j][k];
			}
		}
	}

	for(int i=0; i<nx+2; i++){
		for(int j=0; j<ny+1; j++){
			for(int k=0; k<nz+1; k++){
				VW[i][j][k] = Vaz[i][j][k]*Way[i][j][k];
			}
		}
	}

	return;

}

// consolidate advection terms
void consolidate_advection( boost::multi_array<double, 3>& U,
							boost::multi_array<double, 3>& V,
							boost::multi_array<double, 3>& W,
							boost::multi_array<double, 3>& U2_x,
							boost::multi_array<double, 3>& V2_y,
							boost::multi_array<double, 3>& W2_z,
							boost::multi_array<double, 3>& UV_y,
							boost::multi_array<double, 3>& UW_z,
							boost::multi_array<double, 3>& VU_x,
							boost::multi_array<double, 3>& VW_z,
							boost::multi_array<double, 3>& WU_x,
							boost::multi_array<double, 3>& WV_y,
							cuint nx, cuint ny, cuint nz,
							cdouble dt )
{
	// boost::multi_array<double, 3> U(boost::extents[nx-1][ny][nz]);
	// boost::multi_array<double, 3> V(boost::extents[nx][ny-1][nz]);
	// boost::multi_array<double, 3> W(boost::extents[nx][ny][nz-1]);

	// need to truncate some terms

	// x-direction
	for(int i=0; i<(nx-1); i++){
		for(int j=0; j<(ny); j++){
			for(int k=0; k<(nz); k++){
				U[i][j][k] = U[i][j][k] - dt *
					(U2_x[i][j+1][k+1] + UV_y[i+1][j][k+1] + UW_z[i+1][j+1][k]);
			}
		}
	}

	// y-direction
	for(int i=0; i<(nx); i++){
		for(int j=0; j<(ny-1); j++){
			for(int k=0; k<(nz); k++){
				V[i][j][k] = V[i][j][k] - dt *
					(V2_y[i+1][j][k+1] + VU_x[i][j+1][k+1] + VW_z[i+1][j+1][k]);
			}
		}
	}

	// z-direction
	for(int i=0; i<(nx); i++){
		for(int j=0; j<(ny); j++){
			for(int k=0; k<(nz-1); k++){
				W[i][j][k] = W[i][j][k] - dt *
					(W2_z[i+1][j+1][k] + WU_x[i][j+1][k+1] + WV_y[i+1][j][k+1]);
			}
		}
	}
	
	return;

}
