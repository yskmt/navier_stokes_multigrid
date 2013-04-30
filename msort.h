#include <omp.h>
#include <vector>
#include <iostream>
#include <cstdlib>

using namespace std;

// int(*compar)(const T *, const T *)
// 1 for left<right
// 0 for else

// comparison function for sorting pairs
int comp_tuples( const tuple<uint, uint, double>& i,
				 const tuple<uint, uint, double>& j ) {
    if( (get<0>(i)) < (get<0>(j)) ) return true;
	else if( get<0>(i) == get<0>(j)) return (get<1>(i)) < (get<1>(j));
	else return false;
}

// merge two sorted arrays
void merge(tuple <uint, uint, double>* left,
		   tuple <uint, uint, double>* right,
		   const int n_left, const int n_right,
		   tuple <uint, uint, double>* result,
		   tuple <uint, uint, double>* tmp		   )
{
	unsigned int it = 0;
    unsigned int left_it = 0, right_it = 0;
	// cout<<"n_left "<<n_left<<" n_right "<<n_right<<endl;
	
    while(left_it < n_left && right_it < n_right ) {
		it = left_it+right_it;
		// cout<<it<<endl;
		if(comp_tuples(left[left_it], right[right_it])) {
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

// mergesort with OpenMP parallelism
void mergesort(tuple <uint, uint, double>* vec,
			   const int threads,
			   const int n,
			   tuple <uint, uint, double>* tmp
			   )
{
    // Termination condition: List is completely sorted if it
    // only contains a single element.
    if(n == 1){
		return;
	}

    // Determine the location of the middle element in the vector
	tuple <uint, uint, double>* left = vec; // left array pointer
	int n_left = n/2; // number of elements in left array
	tuple <uint, uint, double>* tmp_left = tmp; // left tmp array pointer
	
	tuple <uint, uint, double>* right = left+n/2; // right array pointer
	int n_right = n-n/2; // number of elements in right array
	tuple <uint, uint, double>* tmp_right = tmp_left+n/2; // right tmp array pointer

    // Perform a merge sort on the two smaller vectors
    if (threads > 1) {
		
		#pragma omp parallel sections
		{
			#pragma omp section
			{
				mergesort(left, threads/2, n_left, tmp_left);
			}
			#pragma omp section
			{
				mergesort(right, threads - threads/2, n_right, tmp_right);
			}
		}
	}
    else {
		mergesort(left, 1, n_left, tmp_left);
		mergesort(right, 1, n_right, tmp_right);
	}

    merge(left, right, n_left, n_right, left, tmp );
	
	return;
}
