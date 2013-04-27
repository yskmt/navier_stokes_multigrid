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

void get_neighbor( unsigned int t[][3][3],
				   cuint i, cuint j, cuint k,
				   cuint I, cuint J, cuint K )
{
	for(int p=0; p<3; p++){
		for(int q=0; q<3; q++){
			for(int r=0; r<3; r++){
				three_d_to_one_d(i-1+p,j-1+q,k-1+r, I,J, t[p][q][r]);
			}
		}
	}
}
