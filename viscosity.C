#include "viscosity.h"


// implicitly solve viscosity
void viscosity(  double* U,
				 double* V,
				 double* W,
				 double* Uss, double* Vss, double* Wss,
				 cuint nx, cuint ny, cuint nz,
				 cdouble hx, cdouble hy, cdouble hz,
				 cdouble hx2i, cdouble hy2i, cdouble hz2i,
				 cdouble dt, cdouble nu,
				 cdouble bcs[][6],
				 cdouble tol, cuint max_iteration )
{
	double Er;
	
	// Lu
	cuint n_u_dof = (nx-1)*ny*nz;
	vector<tuple <uint, uint, double> > Lu_sp;
	vector<double> Lu_val(Lu_sp.size(),0.0);
	vector<uint> Lu_col_ind(Lu_sp.size(), 0);
	vector<uint> Lu_row_ptr(1,0);		
	// double* Fu = new double[n_u_dof];
	// set load vector
	// viscosity_load_vector(Fu, U, nx-1, ny, nz);
	// sparse viscosity matrix and bc modification
	viscosity_matrix_sparse( Lu_sp, Lu_val, Lu_col_ind, Lu_row_ptr,
							 U, nx-1, ny, nz, hx, hy, hz,
							 hx2i, hy2i, hz2i, dt, nu, 
							 bcs[0], X_DIR );
	// now solve Lu\Fu
	// construct solution vector
	double* Uss_tmp = new double[n_u_dof];
	// initial guess
#pragma omp parallel for shared(Uss, Uss_tmp) num_threads(nt)
	for(int n=0; n<n_u_dof; n++){
	    Uss[n] = 0.0;
	    Uss_tmp[n] = 0.0;
    }
	// residual and error
	double* Ru = new double[n_u_dof];
	Er = tol*10;
	// jacobi iteration
	jacobi_sparse(tol, max_iteration, n_u_dof, Uss, Uss_tmp,
				  Lu_val, Lu_col_ind, Lu_row_ptr, U, Er, Ru);
			
	// lv
	cuint n_v_dof = (nx)*(ny-1)*nz;
	vector<tuple <uint, uint, double> > Lv_sp;
	vector<double> Lv_val(Lv_sp.size(),0.0);
	vector<uint> Lv_col_ind(Lv_sp.size(), 0);
	vector<uint> Lv_row_ptr(1,0);		
	// double* Fv = new double[n_v_dof];
	// set load vector
	// viscosity_load_vector(Fv, V, nx, ny-1, nz);
   	// sparse viscosity matrix and bc modification
	viscosity_matrix_sparse( Lv_sp, Lv_val, Lv_col_ind, Lv_row_ptr,
							 V, nx, ny-1, nz, hx, hy, hz,
							 hx2i, hy2i, hz2i, dt, nu, 
							 bcs[1], Y_DIR );
	// now solve Lv\Fv
	// construct solution vector
	// double* Vss = new double[n_v_dof];
	double* Vss_tmp = new double[n_v_dof];
	// initial guess
#pragma omp parallel for shared(Vss, Vss_tmp) num_threads(nt)
	for(int n=0; n<n_v_dof; n++){
	    Vss[n] = 0.0;
	    Vss_tmp[n] = 0.0;
    }
	// residual and error
	double* Rv = new double[n_v_dof];
	Er = tol*10;
	// jacobi iteration
	jacobi_sparse(tol, max_iteration, n_v_dof, Vss, Vss_tmp,
				  Lv_val, Lv_col_ind, Lv_row_ptr, V, Er, Rv);

	// Lw
	cuint n_w_dof = (nx)*(ny)*(nz-1);
	vector<tuple <uint, uint, double> > Lw_sp;
	vector<double> Lw_val(Lw_sp.size(),0.0);
	vector<uint> Lw_col_ind(Lw_sp.size(), 0);
	vector<uint> Lw_row_ptr(1,0);		
	// double* Fw = new double[n_w_dof];
   	// set load vector
	// viscosity_load_vector(Fw, W, nx, ny, nz-1);
	// sparse viscosity matrix and bc modification
	viscosity_matrix_sparse( Lw_sp, Lw_val, Lw_col_ind, Lw_row_ptr,
							 W, nx, ny, nz-1, hx, hy, hz,
							 hx2i, hy2i, hz2i, dt, nu, 
							 bcs[2], Z_DIR );
	// now solve Lw\Fw
	// construct solution vector
	// double* Wss = new double[n_w_dof];
	double* Wss_tmp = new double[n_w_dof];
	// initial guess
#pragma omp parallel for shared(Wss, Wss_tmp) num_threads(nt)
	for(int n=0; n<n_w_dof; n++){
	    Wss[n] = 0.0;
	    Wss_tmp[n] = 0.0;
    }
	// residual and error
	double* Rw = new double[n_w_dof];
	Er = tol*10;
	// jacobi iteration
	jacobi_sparse(tol, max_iteration, n_w_dof, Wss, Wss_tmp,
				  Lw_val, Lw_col_ind, Lw_row_ptr, W, Er, Rw);
	
	// v_cycle( n_u_dof, nx-1, ny, nz, hx2i, hy2i, hz2i,
	// 		 tol, max_iteration, pre_smooth_iteration,
	// 		 lx, ly, lz, 0, max_level, Fu, Er);

	// cleanup
	delete[] Uss_tmp, Vss_tmp, Wss_tmp;
	delete[] Ru, Rv, Rw;
	
}

// sparse viscosity matrix
void viscosity_matrix_sparse( vector<tuple <uint, uint, double> >& L_sp,
							  vector<double>& val,
							  vector<uint>& col_ind,
							  vector<uint>& row_ptr,
							  double* F,
							  cuint nx, cuint ny, cuint nz,
							  cdouble hx, cdouble hy, cdouble hz,
							  cdouble hx2i, cdouble hy2i, cdouble hz2i,
							  cdouble dt, cdouble nu,
							  cdouble* u_bc,
							  cuint dir // direction of flow: u, v, or w?
							  )
{
	// initialize sparse matrix (row#, col#, value)
	vector<vector<tuple <uint, uint, double> > > M;
	M.resize(nt);
	
#pragma omp parallel  shared(M) num_threads(nt)
	{
		cuint myrank = omp_get_thread_num();
		
#pragma omp for 
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			for(int k=0; k<nz; k++){

				unsigned int p,q;
				unsigned int t_011,t_111,t_211,t_101,t_121,t_110,t_112;
				three_d_to_one_d(i,  j,  k,   nx,ny, t_111);
				three_d_to_one_d(i-1,j,  k,   nx,ny, t_011);
				three_d_to_one_d(i+1,j,  k,   nx,ny, t_211);
				three_d_to_one_d(i,  j-1,k,   nx,ny, t_101);
				three_d_to_one_d(i,  j+1,k,   nx,ny, t_121);
				three_d_to_one_d(i,  j,  k-1, nx,ny, t_110);
				three_d_to_one_d(i,  j,  k+1, nx,ny, t_112);
				
				// assignning values
				//U** contribution from left hand side
				sparse_add(M[myrank], t_111, t_111, -1);

				// avoid boundaries
				if(i-1>=0)
					sparse_add(M[myrank], t_111, t_011, hx2i);
				else{ // x0
					if(dir==X_DIR)
						F[t_111] -= dt*nu*u_bc[0]/(hx*hx);
					else{
						F[t_111] -= dt*nu/(hx*hx)*u_bc[0]*2;
						sparse_add(M[myrank], t_111, t_111, -1*hx2i);
					}
				}
					
				sparse_add(M[myrank], t_111, t_111, -2*hx2i);

				if(i+1<nx)
					sparse_add(M[myrank], t_111, t_211, hx2i);
				else{ //xl
					if(dir==X_DIR)
						F[t_111] -= dt*nu/(hx*hx) * u_bc[1];
					else{
						F[t_111] -= dt*nu/(hx*hx)*u_bc[1]*2;
						sparse_add(M[myrank], t_111, t_111, -1*hx2i);
					}
				}
				
				if(j-1>=0)
					sparse_add(M[myrank], t_111, t_101, hy2i);
				else{ // y0
					if(dir==Y_DIR)
						F[t_111] -= dt*nu/(hy*hy) * u_bc[2];
					else{
						F[t_111] -= dt*nu/(hy*hy)*u_bc[2]*2;
						sparse_add(M[myrank], t_111, t_111, -1*hy2i);
					}
						
				}

				sparse_add(M[myrank], t_111, t_111, -2*hy2i);
				
				if(j+1<ny)
					sparse_add(M[myrank], t_111, t_121, hy2i);
				else{ //yl
					if(dir==Y_DIR)
						F[t_111] -= dt*nu/(hy*hy) * u_bc[3];
					else{
						F[t_111] -= dt*nu/(hy*hy)*u_bc[3]*2;
						sparse_add(M[myrank], t_111, t_111, -1*hy2i);
					}
				}
					
				if(k-1>=0)
					sparse_add(M[myrank], t_111, t_110, hz2i);
				else{ // z0
					if(dir==Z_DIR)
						F[t_111] -= dt*nu/(hz*hz) * u_bc[4];
					else{
						F[t_111] -= dt*nu/(hz*hz)*u_bc[4]*2;
						sparse_add(M[myrank], t_111, t_111, -1*hz2i);
					}
				}
					
				sparse_add(M[myrank], t_111, t_111, -2*hz2i);
				
				if(k+1<nz)
					sparse_add(M[myrank], t_111, t_112, hz2i);
				else{ // zl
					if(dir==Z_DIR)
						F[t_111] -= dt*nu/(hz*hz) * u_bc[5];
					else{
						F[t_111] -= dt*nu/(hz*hz)*u_bc[5]*2;
						sparse_add(M[myrank], t_111, t_111, -1*hz2i);
					}
				}
			}
		}
	} // end for

	} // end parallel region		

	// merge and sort
	// cout<<"sorting..."<<endl;
	for(int i=1; i<nt; i++)
		M[0].insert( M[0].end(), M[i].begin(), M[i].end() );
	// sort(M[0].begin(), M[0].end(), comp_pairs);
	vector<tuple <uint, uint, double> > tmp;
	tmp.resize(M[0].size());	
	mergesort(&M[0][0], nt, M[0].size(), &tmp[0] );

	// consolidate
	L_sp.push_back(M[0][0]);
	uint ct=0;
	for(int i =1; i<M[0].size(); i++){
		if( (get<0>(L_sp[ct])==get<0>(M[0][i]))
			&& (get<1>(L_sp[ct])==get<1>(M[0][i])) ){
			get<2>(L_sp[ct]) += get<2>(M[0][i]);
		}
		else{
			L_sp.push_back(M[0][i]);
			ct++;
		}
	}
   
	// convert to CSR format
	// cout<<"converting to CSR format"<<endl;
	val.resize(L_sp.size(),0.0);
	col_ind.resize(L_sp.size(), 0);
	
#pragma omp parallel for shared(val, col_ind, L_sp) num_threads(nt)
	for(int i=0; i<L_sp.size(); i++){
		val[i] = get<2>(L_sp[i]);
		col_ind[i] = get<1>(L_sp[i]);
	}
	for(int i=1; i<L_sp.size(); i++){
		if(get<0>(L_sp[i])!=get<0>(L_sp[i-1]))
		   row_ptr.push_back(i);
	}
	row_ptr.push_back(L_sp.size());

	// cout<<"done"<<endl;
	
	// output to file for testing purpose
	// char file_name[100];
	// if(dir==X_DIR) sprintf(file_name, "L%s_matrix.dat", "u");
	// if(dir==Y_DIR) sprintf(file_name, "L%s_matrix.dat", "v");
	// if(dir==Z_DIR) sprintf(file_name, "L%s_matrix.dat", "w");

	// ofstream file_out(file_name);
	// for(int i=0; i<L_sp.size(); i++){
	// 	file_out<<get<0>(L_sp[i])<<" "<<get<1>(L_sp[i])
	// 		<<" "<<get<2>(L_sp[i])<<endl;
	// }
	// file_out.close();
	
}

// set load vector for implicit viscous solve
void viscosity_load_vector( double* F, double* U,
							cuint nx, cuint ny, cuint nz)
{
	uint t;
#pragma omp parallel for private(t) shared(F, U) num_threads(nt)
	for(int i=0; i<(nx); i++){
		for(int j=0; j<(ny); j++){
			for(int k=0; k<(nz); k++){
				three_d_to_one_d(i,  j,  k, nx,ny, t);

				F[t] = U[t];
			}
		}
	}

	return;
}
