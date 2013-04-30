#include "assemble.h"
#include "utils.h"

// 2nd order stencil
void fd_matrix( double** M,
				cuint I, cuint J, cuint K,
				const double dx2i,
				const double dy2i,
				const double dz2i,
				cuint n_dof				)
{
#pragma omp parallel for shared(M)
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
void fd_matrix_sparse( cuint I, cuint J, cuint K,
					   const double dx2i,
					   const double dy2i,
					   const double dz2i,
					   cuint n_dof )
{
	// initialize sparse matrix (row#, col#, value)
	vector<tuple <uint, uint, double> >M;
	cuint nt=2;

	vector< vector<tuple <uint, uint, double> > > M_sp;
	M_sp.resize(nt);
	
#pragma omp parallel private(M) shared(M_sp) num_threads(nt)
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
				sparse_insert(M, t_111, t_011, dx2i);
				sparse_insert(M, t_111, t_111 , -2*dx2i);
				sparse_insert(M, t_111, t_211 , dx2i);
				// J
				sparse_insert(M, t_111, t_101 , dy2i);
				sparse_insert(M, t_111, t_111 , -2*dy2i);
				sparse_insert(M, t_111, t_121 , dy2i);
				// K
				sparse_insert(M, t_111, t_110 , dz2i);
				sparse_insert(M, t_111, t_111 , -2*dz2i);
				sparse_insert(M, t_111, t_112 , dz2i);

			}
		}
	} // end for

	// sort and consolidate sparse matrix (row#, col#, value)
	// cout<<"sorting..."<<endl;
	sort(M.begin(), M.end(), comp_pairs);
	// cout<<"sorting done"<<endl;
	// vector<tuple <uint, uint, double> >M_sp;
	M_sp[myrank].push_back(M[0]);
	uint ct=0;
		
	for(int i=1; i<M.size(); i++){
		if( (get<0>(M_sp[myrank][ct])==get<0>(M[i]))
			&& (get<1>(M_sp[myrank][ct])==get<1>(M[i])) ){
			// get<0>(M_sp[ct]) += get<0>(M[i]);
			// get<1>(M_sp[ct]) += get<1>(M[i]);
			get<2>(M_sp[myrank][ct]) += get<2>(M[i]);
		}
		else{
			M_sp[myrank].push_back(M[i]);
			ct++;
		}
		// cout<<"i: "<<get<0>(M[i])<<" j: "<<get<1>(M[i])<<" v: "
		// 	<<get<2>(M[i])<<endl;
	}
	
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

	vector<tuple <uint, uint, double> > merged_array, tmp_array;
	merged_array.resize(M_sp[0].size()+M_sp[1].size());
	tmp_array.resize(M_sp[0].size()+M_sp[1].size());
	
	// merge two sorted arrays
   	merge(M_sp[0], M_sp[1],
		  (cuint) M_sp[0].size(), (cuint) M_sp[1].size(),
		  merged_array,
		  tmp_array	
		  );

	// consolidate again...
	vector<tuple <uint, uint, double> >merged_sorted_array;
	merged_sorted_array.push_back(merged_array[0]);
	uint ct=0;
	for(int i=1; i<merged_array.size(); i++){
		if( (get<0>(merged_sorted_array[ct])==get<0>(merged_array[i]))
			&& (get<1>(merged_sorted_array[ct])==get<1>(merged_array[i])) ){
			// get<0>(M_sp[ct]) += get<0>(merged_array[i]);
			// get<1>(M_sp[ct]) += get<1>(merged_array[i]);
			get<2>(merged_sorted_array[ct]) += get<2>(merged_array[i]);
		}
		else{
			merged_sorted_array.push_back(merged_array[i]);
			ct++;
		}
		// cout<<"i: "<<get<0>(M[i])<<" j: "<<get<1>(M[i])<<" v: "
		// 	<<get<2>(M[i])<<endl;
	}
	
	for(int i=0; i<merged_sorted_array.size(); i++){
		cout<<get<0>(merged_sorted_array[i])<<" "<<get<1>(merged_sorted_array[i])
			<<" "<<get<2>(merged_sorted_array[i])<<endl;
	}

	
	// cout<<"setting global constraint"<<endl;
	// global constraint
	// for(int i=0; i<(n_dof); i++){
	// 	M[i][n_dof-1] = 1;
	// 	M[n_dof-1][i] = 1;
	// }
	// M[n_dof-1][n_dof-1] = n_dof;

}


void load_vector( double* F,
				  cuint n_dof,
				  cuint I,
				  cuint J,
				  cuint K
				  )
{
	// construct load vector
	#pragma omp parallel for shared(F)
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
