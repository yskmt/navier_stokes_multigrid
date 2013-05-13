#include "pressure.h"
#include "v_cycle.h"

// compute pressure correction
void pressure( 	double* U, double* V, double* W, double* P,
					double* Uss, double* Vss, double* Wss,
					cuint nx, cuint ny, cuint nz,
					cdouble bcs[][6],
					cdouble lx, cdouble ly, cdouble lz,
					cdouble hx, cdouble hy, cdouble hz,
					cdouble hx2i, cdouble hy2i, cdouble hz2i,
					cdouble tol, cuint max_iteration,
					cuint pre_smooth_iteration, cuint max_level,
					cdouble dt)
{
	cuint n_dof = nx*ny*nz;

	// residual and error
	double* Rp = new double[n_dof];
	double Er = tol*10;
	
	// 0-level v_cycle
	if(max_level==0)
		v_cycle_0( P, Rp,
				   n_dof, nx, ny, nz,
				   hx, hy, hz,
				   hx2i, hy2i, hz2i,
				   tol, max_iteration, pre_smooth_iteration,
				   hx, hy, hz,
				   0, max_level,
				   Er,
				   Uss, Vss, Wss,
				   bcs,dt  );
	else
		// v-cycle
		v_cycle( P, n_dof, nx, ny, nz,
				 hx, hy, hz,
				 hx2i, hy2i, hz2i,
				 tol, max_iteration, pre_smooth_iteration,
				 lx, ly, lz,
				 0, max_level-1,
				 Rp,
				 Er,
				 Uss, Vss, Wss,
				 bcs, dt
				 );

	
	// compute pressure corrections
	double* Pr_x = new double[(nx-1)*(ny)*(nz)];
	double* Pr_y = new double[(nx)*(ny-1)*(nz)];
	double* Pr_z = new double[(nx)*(ny)*(nz-1)];

	for(int i=0; i<(nx-1)*(ny)*(nz); i++)
		Pr_x[i]=0.0;

	compute_corrections(  P, Pr_x, Pr_y, Pr_z, nx, ny, nz, hy, hy, hz );
	
	// 1d index
	uint t;

	// correct velocities
	// x-direction
#pragma omp parallel for private(t) shared(U, Uss, Pr_x) num_threads(nt)
	for(int i=0; i<nx-1; i++){
		for(int j=0; j<ny; j++){
			for(int k=0; k<nz; k++){
				three_d_to_one_d(i,j,k, nx-1, ny, t);
				U[t] = Uss[t] - Pr_x[t]*dt;
			}
		}
	}
	
	// y-direction	
#pragma omp parallel for private(t) shared(V, Vss, Pr_y) num_threads(nt)
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny-1; j++){
			for(int k=0; k<nz; k++){
				three_d_to_one_d(i,j,k, nx, ny-1, t);
				V[t] = Vss[t] - Pr_y[t]*dt;
			}
		}
	}
	
	// z-direction
#pragma omp parallel for private(t) shared(W, Wss, Pr_z) num_threads(nt)
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			for(int k=0; k<nz-1; k++){
				three_d_to_one_d(i,j,k, nx, ny, t);
				W[t] = Wss[t] - Pr_z[t]*dt;
			}
		}
	}

	// cleanup
	delete[] Pr_x, Pr_y, Pr_z;

}

// build right hand side of pressure poisson equation
void pressure_rhs( double* F,
				   double* Uss, double* Vss, double* Wss,
				   cuint nx, cuint ny, cuint nz,
				   cdouble bcs[][6],
				   cdouble hx, cdouble hy, cdouble hz,
				   cdouble dt
				   )
{
	// for(int i=0; i<(nx-1)*ny*nz; i++)
	// 	cout<<"Uss: "<<Uss[i]<<endl;
	// for(int i=0; i<(nx)*(ny-1)*nz; i++)
	// 	cout<<"Vss: "<<Vss[i]<<endl;
	// for(int i=0; i<(nx)*ny*(nz-1); i++)
	// 	cout<<"Wss: "<<Wss[i]<<endl;
	cuint n_dof = nx*ny*nz;

	// initialize
	for(int i=0; i<n_dof; i++)
		F[i] =0;
	
	// Uss contribution
#pragma omp parallel for shared(Uss, F) num_threads(nt)
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			for(int k=0; k<nz; k++){
				uint t0, t1, t;
				three_d_to_one_d(i-1,j,k, nx-1,ny, t0);
				three_d_to_one_d(i,j,k, nx-1,ny, t1);
				three_d_to_one_d(i,j,k, nx,ny, t);
				if(i==0)
					F[t] += (Uss[t1]-bcs[0][0])/hx;
				else if(i==nx-1)
					F[t] += (bcs[0][1]-Uss[t0])/hx;
				else{
					F[t] += (Uss[t1]-Uss[t0]) / hx;
				}
			}
		}
	}

	// Vss contribution
#pragma omp parallel for shared(Vss, F) num_threads(nt)
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			for(int k=0; k<nz; k++){
				uint t0, t1, t;
				three_d_to_one_d(i,j-1,k, nx,ny-1, t0);
				three_d_to_one_d(i,j,k, nx,ny-1, t1);
				three_d_to_one_d(i,j,k, nx,ny, t);
				if(j==0)
					F[t] += (Vss[t1]-bcs[1][2])/hy;
				else if(j==ny-1)
					F[t] += (bcs[1][3]-Vss[t0])/hy;
				else{
					F[t] += (Vss[t1]-Vss[t0]) / hy;
				}
			}
		}
	}
	
	// Wss contribution
#pragma omp parallel for shared(Wss, F) num_threads(nt)
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			for(int k=0; k<nz; k++){
				uint t0, t1, t;
				three_d_to_one_d(i,j,k-1, nx,ny, t0);
				three_d_to_one_d(i,j,k, nx,ny, t1);
				three_d_to_one_d(i,j,k, nx,ny, t);
				if(k==0)
					F[t] += (Wss[t1]-bcs[2][4])/hz;
				else if(k==nz-1){
					F[t] += (bcs[2][5]-Wss[t0])/hz;
				}
				else{
					F[t] += (Wss[t1]-Wss[t0]) / hz;
				}
			}
		}
	}

	// divide everyhing by dt
	for(int i=0; i<(nx*ny*nz); i++){
		F[i] = F[i]/dt;
	}
									 
	
	// global constraint to close the system
	// F[n_dof] = 0;

	// point boundary condition at the center of domain
	// uint t;
	// three_d_to_one_d(uint(nx/2),uint(ny/2),uint(nz/2), nx,ny, t);
	// F[t] = 0;
	
	// output to file for testing purpose
	// ofstream file_out("Fp_vector.dat");
	// for(int i=0; i<(n_dof); i++){
	// 	file_out<<F[i]<<endl;
	// }
	// file_out.close();

	
	return;
}

// 2nd order stencil
void pressure_matrix( vector<tuple <uint, uint, double> >& Lp_sp,
					  vector<double>& val,
					  vector<uint>& col_ind,
					  vector<uint>& row_ptr,
					  cuint nx, cuint ny, cuint nz,
					  const double hx2i,
					  const double hy2i,
					  const double hz2i,
					  cuint n_dof
					  )
{
	// initialize sparse matrix (row#, col#, value)
	vector<vector<tuple <uint, uint, double> > > M;
	M.resize(nt);
	
#pragma omp parallel  shared(M) num_threads(nt)
	{
		cuint myrank = omp_get_thread_num();

		// loop through inner nodes
#pragma omp for 
		for(int i=0; i<nx; i++){
			for(int j=0; j<ny; j++){
				for(int k=0; k<nz; k++){
					// 2nd order stencil
					unsigned int p,q;
					unsigned int t_011,t_111,t_211,t_101,t_121,t_110,t_112;
					three_d_to_one_d(i,  j,  k,   nx,ny, t_111);
					three_d_to_one_d(i-1,j,  k,   nx,ny, t_011);
					three_d_to_one_d(i+1,j,  k,   nx,ny, t_211);
					three_d_to_one_d(i,  j-1,k,   nx,ny, t_101);
					three_d_to_one_d(i,  j+1,k,   nx,ny, t_121);
					three_d_to_one_d(i,  j,  k-1, nx,ny, t_110);
					three_d_to_one_d(i,  j,  k+1, nx,ny, t_112);

					// p_xx contributions
					if(i-1>=0)
						sparse_add(M[myrank], t_111, t_011, hx2i);
					else // x0: i==0 (P[-1][j][k]==P[0][j][k])
						sparse_add(M[myrank], t_111, t_111, hx2i);
					
					sparse_add(M[myrank], t_111, t_111, -2*hx2i);

					if(i+1<nx)
						sparse_add(M[myrank], t_111, t_211, hx2i);
					else // xl: i==nx-1
						sparse_add(M[myrank], t_111, t_111, hx2i);

					// p_yy contributions
					if(j-1>=0)
						sparse_add(M[myrank], t_111, t_101, hy2i);
					else // y0: j==0 (P[i][-1][k]==P[i][0][k])
						sparse_add(M[myrank], t_111, t_111, hy2i);
					
					sparse_add(M[myrank], t_111, t_111, -2*hy2i);

					if(j+1<ny)
						sparse_add(M[myrank], t_111, t_121, hy2i);
					else // yl: j==ny-1
						sparse_add(M[myrank], t_111, t_111, hy2i);

					// p_zz contributions
					if(k-1>=0)
						sparse_add(M[myrank], t_111, t_110, hz2i);
					else // z0: k==0 (P[i][j][-1]==P[i][j][0])
						sparse_add(M[myrank], t_111, t_111, hz2i);
					
					sparse_add(M[myrank], t_111, t_111, -2*hz2i);

					if(k+1<nz)
						sparse_add(M[myrank], t_111, t_112, hz2i);
					else // zl: k==nz-1
						sparse_add(M[myrank], t_111, t_111, hz2i);

				}
			}
		} // end for


		// global constraint to close the system
		// #pragma omp for
		// 		for(int i=0; i<n_dof; i++){
		// 			sparse_add(M[myrank], i, n_dof, 1);
		// 			sparse_add(M[myrank], n_dof, i, 1);
		// 		}
		// 		if(myrank==0)
		// 			sparse_add(M[myrank], n_dof, n_dof,
		// 					   n_dof);
		
	} // end parallel region		

	// merge and sort
	for(int i=1; i<nt; i++)
		M[0].insert( M[0].end(), M[i].begin(), M[i].end() );
	vector<tuple <uint, uint, double> > tmp;
	tmp.resize(M[0].size());	
	mergesort(&M[0][0], nt, M[0].size(), &tmp[0] );
	
	// consolidate
	Lp_sp.push_back(M[0][0]);
	uint ct=0;
	for(int i =1; i<M[0].size(); i++){
		if( (get<0>(Lp_sp[ct])==get<0>(M[0][i]))
			&& (get<1>(Lp_sp[ct])==get<1>(M[0][i])) ){
			get<2>(Lp_sp[ct]) += get<2>(M[0][i]);
		}
		else{
			Lp_sp.push_back(M[0][i]);
			ct++;
		}
	}

	// point constraint to close the system
	get<2>(Lp_sp[0]) = 5*get<2>(Lp_sp[0]);
		
	// convert to CSR format
	val.resize(Lp_sp.size(),0.0);
	col_ind.resize(Lp_sp.size(), 0);
	
#pragma omp parallel for shared(val, col_ind, Lp_sp) num_threads(nt)
	for(int i=0; i<Lp_sp.size(); i++){
		val[i] = get<2>(Lp_sp[i]);
		col_ind[i] = get<1>(Lp_sp[i]);
	}
	for(int i=1; i<Lp_sp.size(); i++){
		if(get<0>(Lp_sp[i])!=get<0>(Lp_sp[i-1]))
			row_ptr.push_back(i);
	}
	row_ptr.push_back(Lp_sp.size());

	// output to file for testing purpose
	ofstream file_out("Lp_matrix.dat");
	for(int i=0; i<Lp_sp.size(); i++){
		file_out<<get<0>(Lp_sp[i])<<" "<<get<1>(Lp_sp[i])
				<<" "<<get<2>(Lp_sp[i])<<endl;
	}
	file_out.close();
	
}

// compute corrections from pressure value
void compute_corrections( double* Pr,
						  double* Pr_x,
						  double* Pr_y,
						  double* Pr_z,
						  cuint nx, cuint ny, cuint nz,
						  cdouble hx, cdouble hy, cdouble hz )
{
	for(int i=0; i<(nx-1)*(ny)*(nz); i++)
		Pr_x[i]=0.0;
	for(int i=0; i<(nx)*(ny-1)*(nz); i++)
		Pr_y[i]=0.0;
	for(int i=0; i<(nx)*(ny)*(nz-1); i++)
		Pr_z[i]=0.0;
	
	staggered_first_difference( Pr, Pr_x, nx, ny, nz, nx-1, ny, nz, hx, X_DIR );
	staggered_first_difference( Pr, Pr_y, nx, ny, nz, nx, ny-1, nz, hy, Y_DIR );
	staggered_first_difference( Pr, Pr_z, nx, ny, nz, nx, ny, nz-1, hz, Z_DIR );

	return;
}
	
