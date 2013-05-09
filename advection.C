#include "advection.h"

// set up initial conditions
void initial_conditoins()
{

	return;
	
}

// treat nonlinear advection terms
void advection( cuint nx, cuint ny, cuint nz )
{
	// boundary values
	double uB = 0.0; //bottom
	double uT = 0.0; //top
	
	boost::multi_array<double, 3> Ue(boost::extents[nx+1][ny+2][nz+2]);
	boost::multi_array<double, 3> Ve(boost::extents[nx+2][ny+1][nz+2]);
	boost::multi_array<double, 3> We(boost::extents[nx+2][ny+2][nz+1]);

	grid_matrix(Ue, nx, ny, nz, X_DIR, 0,0,0,0, uB, uT);
	grid_matrix(Ve, nx, ny, nz, Y_DIR, 0,0,0,0, uB, uT);
	grid_matrix(We, nx, ny, nz, Z_DIR, 0,0,0,0, uB, uT);
	
	// get (U^2),x defined at cell centers
	boost::multi_array<double, 3> Ua(boost::extents[nx+1][ny+1][nz+1]);
	boost::multi_array<double, 3> Va(boost::extents[nx+1][ny+1][nz+1]);
	boost::multi_array<double, 3> Wa(boost::extents[nx+1][ny+1][nz+1]);

	// average into cell vertices
	average(Ue, Ua, YZ_DIR);
	average(Ve, Va, XZ_DIR);
	average(We, Wa, XY_DIR);
	

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


}
