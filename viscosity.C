#include "viscosity.h"


// implicitly solve viscosity
void viscosity( cuint nx, cuint ny, cuint nz,
				cdouble hx2i, cdouble hy2i, cdouble hz2i,
				cdouble dt, cdouble nu )
{
	// Lu
	cuint n_u_dof = (nx-1)*ny*nz;
	vector<tuple <uint, uint, double> > Lu_sp;
	vector<double> Lu_val(Lu_sp.size(),0.0);
	vector<uint> Lu_col_ind(Lu_sp.size(), 0);
	vector<uint> Lu_row_ptr(1,0);		
	viscosity_matrix_sparse( Lu_sp,
							 Lu_val, Lu_col_ind, Lu_row_ptr,
							 nx-1, ny, nz,
							 hx2i*dt*nu, hy2i*dt*nu, hz2i*dt*nu,
							 n_u_dof,
							 X_DIR);

	// // Lv
	// cuint n_v_dof = (nx)*(ny-1)*nz;
	// vector<tuple <uint, uint, double> > Lv_sp;
	// vector<double> Lv_val(Lv_sp.size(),0.0);
	// vector<uint> Lv_col_ind(Lv_sp.size(), 0);
	// vector<uint> Lv_row_ptr(1,0);		
	// viscosity_matrix_sparse( Lv_sp,
	// 						 Lv_val, Lv_col_ind, Lv_row_ptr,
	// 						 nx, ny, nz,
	// 						 hx2i*dt*nu, hy2i*dt*nu, hz2i*dt*nu,
	// 						 n_v_dof,
	// 						 Y_DIR);
	
	// // Lw
	// cuint n_w_dof = (nx)*ny*(nz-1);
	// vector<tuple <uint, uint, double> > Lw_sp;
	// vector<double> Lw_val(Lw_sp.size(),0.0);
	// vector<uint> Lw_col_ind(Lw_sp.size(), 0);
	// vector<uint> Lw_row_ptr(1,0);		
	// viscosity_matrix_sparse( Lw_sp,
	// 						 Lw_val, Lw_col_ind, Lw_row_ptr,
	// 						 nx, ny, nz,
	// 						 hx2i*dt*nu, hy2i*dt*nu, hz2i*dt*nu,
	// 						 n_w_dof,
	// 						 Z_DIR);

	
}


// sparse viscosity matrix
void viscosity_matrix_sparse( vector<tuple <uint, uint, double> >& L_sp,
							  vector<double>& val,
							  vector<uint>& col_ind,
							  vector<uint>& row_ptr,
							  cuint nx, cuint ny, cuint nz,
							  cdouble hx2i,
							  cdouble hy2i,
							  cdouble hz2i,
							  cuint n_dof,
							  cuint dir // direction of flow: u, v, or w?
							  )
{
	// initialize sparse matrix (row#, col#, value)
	vector<vector<tuple <uint, uint, double> > > M;
	M.resize(nt);

	uint nx_ = nx;
	uint ny_ = ny;
	uint nz_ = nz;
	
	if(dir==X_DIR) nx_ = nx-1;
	else if(dir==Y_DIR) ny_ = ny-1;
	else if(dir==Z_DIR) nz_ = nz-1;
	
#pragma omp parallel  shared(M) num_threads(nt)
	{
		cuint myrank = omp_get_thread_num();
		
#pragma omp for 
	for(int i=1; i<nx_-1; i++){
		for(int j=1; j<ny_-1; j++){
			for(int k=1; k<nz_-1; k++){

				unsigned int p,q;
				unsigned int t_011,t_111,t_211,t_101,t_121,t_110,t_112;
				three_d_to_one_d(i,  j,  k,   nx_,ny_, t_111);
				three_d_to_one_d(i-1,j,  k,   nx_,ny_, t_011);
				three_d_to_one_d(i+1,j,  k,   nx_,ny_, t_211);
				three_d_to_one_d(i,  j-1,k,   nx_,ny_, t_101);
				three_d_to_one_d(i,  j+1,k,   nx_,ny_, t_121);
				three_d_to_one_d(i,  j,  k-1, nx_,ny_, t_110);
				three_d_to_one_d(i,  j,  k+1, nx_,ny_, t_112);
				
				// assignning values
				//U** contribution from left hand side
				sparse_add(M[myrank], t_111, t_111, -1);

				// avoid boundaries
				if(i-1!=0)
					sparse_add(M[myrank], t_111, t_011, hx2i);
				sparse_add(M[myrank], t_111, t_111, -2*hx2i);
				if(i+1!=nx_-1)
					sparse_add(M[myrank], t_111, t_211, hx2i);
				if(j-1!=0)
					sparse_add(M[myrank], t_111, t_101, hy2i);
				sparse_add(M[myrank], t_111, t_111, -2*hy2i);
				if(j+1!=ny_-1)
					sparse_add(M[myrank], t_111, t_121, hy2i);
				if(k-1!=0)
					sparse_add(M[myrank], t_111, t_110, hz2i);
				sparse_add(M[myrank], t_111, t_111, -2*hz2i);
				if(k+1!=nz_-1)
					sparse_add(M[myrank], t_111, t_112, hz2i);

			}
		}
	} // end for

	} // end parallel region		

	// merge and sort
	cout<<"sorting..."<<endl;
	for(int i=1; i<nt; i++)
		M[0].insert( M[0].end(), M[i].begin(), M[i].end() );
	// sort(M[0].begin(), M[0].end(), comp_pairs);
	vector<tuple <uint, uint, double> > tmp;
	tmp.resize(M[0].size());	
	mergesort(&M[0][0], nt, M[0].size(), &tmp[0] );

	cout<<"sorting done"<<endl;
	
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
	cout<<"converting to CSR format"<<endl;
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

	cout<<"done"<<endl;
	
	// output to file for testing purpose
	char file_name[100];
	if(dir==X_DIR) sprintf(file_name, "L%s_matrix.dat", "u");
	if(dir==Y_DIR) sprintf(file_name, "L%s_matrix.dat", "v");
	if(dir==Z_DIR) sprintf(file_name, "L%s_matrix.dat", "w");

	ofstream file_out(file_name);
	for(int i=0; i<L_sp.size(); i++){
		file_out<<get<0>(L_sp[i])<<" "<<get<1>(L_sp[i])
			<<" "<<get<2>(L_sp[i])<<endl;
	}
	file_out.close();
	
}

// set dirichlet boudnary condition for implicit viscous solve
void set_viscouse_bc( const unsigned int n_dof,
					  const unsigned int nx,
					  const unsigned int ny,
					  const unsigned int nz,
					  vector<vector<tuple <uint, uint, double> > >& M,
					  double* F,
					  cuint myrank
					  )
{
	// x0, xl
	for(int j=0; j<ny; j++){
		for(int k=0; k<nz; k++){
			uint t0, tl;
			three_d_to_one_d(0,  j,  k, nx,ny, t0);
			three_d_to_one_d(nx-1,  j,  k, nx,ny, tl);

			sparse_add(M[myrank], t0, t0, 1);
			sparse_add(M[myrank], tl, tl, 1);
		}
	}

	// y0, yl
	for(int i=0; i<nx; i++){
		for(int k=0; k<nz; k++){
			uint t0, tl;
			three_d_to_one_d(i,  0,  k, nx,ny, t0);
			three_d_to_one_d(i,  ny-1,  k, nx,ny, tl);

			sparse_add(M[myrank], t0, t0, 1);
			sparse_add(M[myrank], tl, tl, 1);
		}
	}
	
	// z0, zl
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			uint t0, tl;
			three_d_to_one_d(i,  j,  0, nx,ny, t0);
			three_d_to_one_d(i,  j,  nz-1, nx,ny, tl);

			sparse_add(M[myrank], t0, t0, 1);
			sparse_add(M[myrank], tl, tl, 1);
		}
	}
		
	return;
}
