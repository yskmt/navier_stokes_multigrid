#include "assemble.h"
#include "utils.h"
#include "msort.h"

// 2nd order stencil
void fd_matrix( double** M,
				cuint I, cuint J, cuint K,
				const double dx2i,
				const double dy2i,
				const double dz2i,
				cuint n_dof
				 )
{
#pragma omp parallel for shared(M) num_threads(nt)
	for(int i=0; i<I; i++){
		for(int j=0; j<J; j++){
			for(int k=0; k<K; k++){
				unsigned int p,q;
				unsigned int t_011,t_111,t_211,t_101,t_121,t_110,t_112;
				three_d_to_one_d(i,j,k, I,J, t_111);
				if(i==0)
					three_d_to_one_d(I-1,j,k, I,J, t_011);
				else
					three_d_to_one_d(i-1,j,k, I,J, t_011);
				if(i==(I-1))
					three_d_to_one_d(0,j,k, I,J, t_211);
				else
					three_d_to_one_d(i+1,j,k, I,J, t_211);

				if(j==0)
					three_d_to_one_d(i,J-1,k, I,J, t_101);
				else
					three_d_to_one_d(i,j-1,k, I,J, t_101);
				if(j==(J-1))
					three_d_to_one_d(i,0,k, I,J, t_121);
				else
					three_d_to_one_d(i,j+1,k, I,J, t_121);
								
				if(k==0)
					three_d_to_one_d(i,j,K-1, I,J, t_110);
				else
					three_d_to_one_d(i,j,k-1, I,J, t_110);
				if(k==(K-1))
					three_d_to_one_d(i,j,0, I,J, t_112);
				else
					three_d_to_one_d(i,j,k+1, I,J, t_112);
				
				// I
				M[t_111][t_011] += dx2i;
				M[t_111][t_111] += -2*dx2i;
				M[t_111][t_211] += dx2i;
				// J
				M[t_111][t_101] += dy2i;
				M[t_111][t_111] += -2*dy2i;
				M[t_111][t_121] += dy2i;
				// K
				M[t_111][t_110] += dz2i;
				M[t_111][t_111] += -2*dz2i;
				M[t_111][t_112] += dz2i;

			}
		}
	}

	// cout<<"setting global constraint"<<endl;
	// global constraint
	for(int i=0; i<(n_dof); i++){
		M[i][n_dof-1] = 1;
		M[n_dof-1][i] = 1;
	}
	M[n_dof-1][n_dof-1] = n_dof;

}

// 2nd order stencil
void fd_matrix_sparse( 	vector<tuple <uint, uint, double> >& M_sp,
						vector<double>& val,
						vector<uint>& col_ind,
						vector<uint>& row_ptr,
						cuint I, cuint J, cuint K,
						const double dx2i,
						const double dy2i,
						const double dz2i,
						cuint n_dof
						)
{

	// initialize sparse matrix (row#, col#, value)
	vector<vector<tuple <uint, uint, double> > > M;
	M.resize(nt);

	
#pragma omp parallel  shared(M) num_threads(nt)
	{
		cuint myrank = omp_get_thread_num();
		
#pragma omp for 
	for(int i=0; i<I; i++){
		for(int j=0; j<J; j++){
			for(int k=0; k<K; k++){
				unsigned int p,q;
				unsigned int t_011,t_111,t_211,t_101,t_121,t_110,t_112;
				three_d_to_one_d(i,j,k, I,J, t_111);
				if(i==0)
					three_d_to_one_d(I-1,j,k, I,J, t_011);
				else
					three_d_to_one_d(i-1,j,k, I,J, t_011);
				if(i==(I-1))
					three_d_to_one_d(0,j,k, I,J, t_211);
				else
					three_d_to_one_d(i+1,j,k, I,J, t_211);

				if(j==0)
					three_d_to_one_d(i,J-1,k, I,J, t_101);
				else
					three_d_to_one_d(i,j-1,k, I,J, t_101);
				if(j==(J-1))
					three_d_to_one_d(i,0,k, I,J, t_121);
				else
					three_d_to_one_d(i,j+1,k, I,J, t_121);
								
				if(k==0)
					three_d_to_one_d(i,j,K-1, I,J, t_110);
				else
					three_d_to_one_d(i,j,k-1, I,J, t_110);
				if(k==(K-1))
					three_d_to_one_d(i,j,0, I,J, t_112);
				else
					three_d_to_one_d(i,j,k+1, I,J, t_112);

				// I
				sparse_insert(M[myrank], t_111, t_011, dx2i);
				sparse_insert(M[myrank], t_111, t_111 , -2*dx2i);
				sparse_insert(M[myrank], t_111, t_211 , dx2i);
				// J
				sparse_insert(M[myrank], t_111, t_101 , dy2i);
				sparse_insert(M[myrank], t_111, t_111 , -2*dy2i);
				sparse_insert(M[myrank], t_111, t_121 , dy2i);
				// K
				sparse_insert(M[myrank], t_111, t_110 , dz2i);
				sparse_insert(M[myrank], t_111, t_111 , -2*dz2i);
				sparse_insert(M[myrank], t_111, t_112 , dz2i);

			}
		}
	} // end for


	// cout<<"setting global constraint"<<endl;
	// global constraint
#pragma omp for
	for(int i=0; i<(n_dof-1); i++){
		sparse_insert(M[myrank], i, n_dof-1, 1);
		sparse_insert(M[myrank], n_dof-1, i, 1);
			// M[i][n_dof-1] = 1;
		// M[n_dof-1][i] = 1;
	}
	if(myrank==0)
		sparse_insert(M[myrank], n_dof-1, n_dof-1,
					  n_dof);
	// else
	// 	sparse_insert(M[myrank], n_dof-1, n_dof-1,
	// 				  -get<2>(M[myrank][n_dof-1]));
	// M[n_dof-1][n_dof-1] = n_dof;
	
	// sort and consolidate sparse matrix (row#, col#, value)
	// cout<<"sorting..."<<endl;
	// sort(M[myrank].begin(), M[myrank].end(), comp_pairs);
	// cout<<"sorting done"<<endl;
	// vector<tuple <uint, uint, double> >M_sp;
	// M_sp[myrank].push_back(M[0]);
	uint ct=0;
			
// #pragma omp critical
// 	{
// 		cout<<"thread #: "<<omp_get_thread_num()<<endl;
// 		for(int i=0; i<M_sp[myrank].size(); i++){
// 			cout<<"i: "<<get<0>(M_sp[myrank][i])<<" j: "
// 				<<get<1>(M_sp[myrank][i])<<" v: "
// 				<<get<2>(M_sp[myrank][i])<<endl;
// 		}
// 	}
	
	} // end parallel region		

	// merge and sort
	cout<<"sorting..."<<endl;
	for(int i=1; i<nt; i++)
		M[0].insert( M[0].end(), M[i].begin(), M[i].end() );
	// sort(M[0].begin(), M[0].end(), comp_pairs);
	vector<tuple <uint, uint, double> > tmp;
	tmp.resize(M[0].size());	
	mergesort(&M[0][0], nt, M[0].size(), &tmp[0] );


	cout<<"done"<<endl;
	
	// consolidate
	M_sp.push_back(M[0][0]);
	uint ct=0;
	for(int i =1; i<M[0].size(); i++){
		if( (get<0>(M_sp[ct])==get<0>(M[0][i]))
			&& (get<1>(M_sp[ct])==get<1>(M[0][i])) ){
			// get<0>(M_sp[ct]) += get<0>(M[0][i]);
			// get<1>(M_sp[ct]) += get<1>(M[0][i]);
			get<2>(M_sp[ct]) += get<2>(M[0][i]);
		}
		else{
			M_sp.push_back(M[0][i]);
			ct++;
		}

	}

   
	// convert to CSR format
	cout<<"converting to CSR format"<<endl;
	val.resize(M_sp.size(),0.0);
	col_ind.resize(M_sp.size(), 0);
	
#pragma omp parallel for shared(val, col_ind, M_sp) num_threads(nt)
	for(int i=0; i<M_sp.size(); i++){
		val[i] = get<2>(M_sp[i]);
		col_ind[i] = get<1>(M_sp[i]);
	}
	for(int i=1; i<M_sp.size(); i++){
		if(get<0>(M_sp[i])!=get<0>(M_sp[i-1]))
		   row_ptr.push_back(i);
	}
	row_ptr.push_back(M_sp.size());

	// for(int i=0; i<row_ptr.size(); i++)
	// 	cout<<row_ptr[i]<<endl;
	cout<<"done"<<endl;
	
	// // output to file for testing purpose
	// ofstream file_out("test_sp_matrix.dat");
	// for(int i=0; i<M_sp.size(); i++){
	// 	file_out<<get<0>(M_sp[i])<<" "<<get<1>(M_sp[i])
	// 		<<" "<<get<2>(M_sp[i])<<endl;
	// }
	// file_out.close();
	
}

// assemble a load vector (only on level 0)
void load_vector( double* F,
				  cuint n_dof,
				  cuint I,
				  cuint J,
				  cuint K
				  )
{
	// construct load vector
#pragma omp parallel for shared(F) num_threads(nt)
	for(int n=0; n<n_dof-1; n++){
		unsigned int i,j,k;
		one_d_to_three_d( n, I, J, i, j, k);
	    F[n] = sin(double(i)/double(I)*2*pi); // solution=(2pi)^2*sin(2pi*x);

	    // F[n] = sin(double(i)/double(I)*2*pi) * sin(double(j)/double(J)*2*pi)
			// * sin(double(k)/double(K)*2*pi);
     }

	// global constraint
	F[n_dof-1] = 0;

}

// set dirichlet boudnary condition
// for periodic domain, it is sufficient to set only one point
int boundary_conditins( const unsigned int n_dof,
		const unsigned int I,
		const unsigned int J,
		const unsigned int K,
		double** M,
		double* F
		)
 {
	int n_bd=0;
	// boundary conditions
	uint t;
	three_d_to_one_d(0,0,0, I,J, t);		
	for(int n=0; n<n_dof; n++){
	    M[n][t]=0;
		M[t][n]=0;
	}
	M[t][t] = 1;
	F[t] = 1;
	
	// #pragma omp parallel for shared(M)
	// for(int i=0; i<I; i++){
	// 	for(int j=0; j<J; j++){
	// 		for(int k=0; k<K; k++){
	// 			if(i==0 || j==0 || k==0
	// 				|| i==(I-1) || j==(J-1) || k==(K-1) ){
	// 				n_bd++;
	// 				unsigned int t;
	// 				three_d_to_one_d(i,j,k, I,J, t);
					
	// 				for(int n=0; n<n_dof; n++){
	// 					M[n][t]=0;
	// 					M[t][n]=0;
	// 				}
	// 				M[t][t] = 1;
	// 			}
	// 		}
	// 	}
	// }

	return n_bd;
}
 
// insert index and value into a sparse matrix
 void sparse_insert( vector<tuple<uint, uint, double > >& M,
		cuint i, cuint j, cdouble v)
{
	tuple<uint, uint, double > M_tmp(i, j, v);
	M.push_back(M_tmp);

	// vector<int> idx_tmp(2, 0.0);
	// idx_tmp[0] = i; idx_tmp[1]=j;
	// idx.push_back(idx_tmp);
	// value.push_back(v);
					
}
				
// merge two sorted arrays 
 void merge(vector<tuple <uint, uint, double> >& left,
		vector<tuple <uint, uint, double> >& right,
		cuint n_left, cuint n_right,
		vector<tuple <uint, uint, double> >& result,
		vector<tuple <uint, uint, double> >& tmp
		)
{
	uint it = 0;
    uint left_it = 0, right_it = 0;
	
    while(left_it < n_left && right_it < n_right ) {
		it = left_it+right_it;
		// cout<<it<<endl;
		if(comp_pairs(left[left_it], right[right_it])) {
			tmp[it] = left[left_it];
			left_it++;
		}
		else{
			tmp[it] = right[right_it];
			right_it++;
		}
	}

    // Push the remaining data from both vectors onto the tmp
    while(left_it < n_left) {
		it = left_it+right_it;
        tmp[it] = left[left_it];
        left_it++;
    }

    while(right_it < n_right) {
		it = left_it+right_it;
        tmp[it] = right[right_it];
        right_it++;
    }

	// Finally put everyhing in result array
	for(int i=0; i<(n_right+n_left); i++)
		result[i] = tmp[i];

}
