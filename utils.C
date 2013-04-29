#include "utils.h"

void three_d_to_one_d( const unsigned int i,
					  const unsigned int j,
					  const unsigned int k,
					  const unsigned int I,
					  const unsigned int J,
					  unsigned int& t )
{
	t=i + j*I + k*I*J;
}

void one_d_to_three_d( const unsigned int t,
					   const unsigned int I,
					   const unsigned int J,
					   unsigned int& i,
					   unsigned int& j,
					   unsigned int& k)
{
	k = t/(I*J);
	j = (t-k*I*J)/I;
	i = t-j*I - k*I*J;
}

// get the neighboring node numbers (periodic domain)
// watch out for the negative unsigned int!!
void get_neighbor( uint t[][3][3],
				   cuint i, cuint j, cuint k,
				   cuint I, cuint J, cuint K )
{
	for(int p=0; p<3; p++){
		for(int q=0; q<3; q++){
			for(int r=0; r<3; r++){
				int nei_i, nei_j, nei_k;
				if(i+p>=I+1) nei_i = i+p-1-I;
				else if(i+p<1) nei_i = i+p+I-1; 
				else nei_i = i+p-1;

				if(j+q>=J+1) nei_j = j+q-1-J;
				else if(j+q<1) nei_j = j+q+J-1; 
				else nei_j = j+q-1;
				
				if(k+r>=K+1) nei_k = k+r-1-K;
				else if(k+r<1) nei_k = k+r+K-1; 
				else nei_k = k+r-1;
				// cout<<nei_i<<" "<<nei_j<<" "<<nei_k<<endl;
				
				// int nei_i = ((i+p)<1) ? (i+p+I-1) : (i+p-1);
				// int nei_j = ((j+p)<1) ? (j+p+J-1) : (j+p-1);
				// int nei_k = ((k+p)<1) ? (k+p+K-1) : (k+p-1);
				three_d_to_one_d(nei_i,nei_j,nei_k, I,J, t[p][q][r]);
			}
		}
	}
}

// get the node numbers in a box (periodic domain)
// watch out for the negative unsigned int!!
void get_box( uint t[][2][2],
				   cuint i, cuint j, cuint k,
				   cuint I, cuint J, cuint K )
{
	for(int p=0; p<2; p++){
		for(int q=0; q<2; q++){
			for(int r=0; r<2; r++){
				int nei_i, nei_j, nei_k;
				if(i+p>=I) nei_i = i+p-I;
				else nei_i = i+p;
				if(j+q>=J) nei_j = j+q-J;
				else nei_j = j+q;
				if(k+r>=K) nei_k = k+r-K;
				else nei_k = k+r;
				
				three_d_to_one_d(nei_i,nei_j,nei_k, I,J, t[p][q][r]);
			}
		}
	}
}
