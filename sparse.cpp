#include "mmreader.hpp"
#include <time.h>
#include <cstring>
#include <iostream>
#include <sys/time.h>
#include <unistd.h>
#include <algorithm>
#include <limits>

using namespace std;

struct args {
	struct sparse_mtx *A;
	struct dense_mtx *B;
	struct dense_mtx *C;
	int p;
	int t_num;
};
pthread_mutex_t C_lock = PTHREAD_MUTEX_INITIALIZER;

void *mul_sparse_dense(void *arg);
int find_row(struct sparse_mtx *A, int val_i);
bool nearly_equal(float a, float b, float epsilon);


bool
SCsrMatrixfromFile(struct sparse_mtx *A, const char* filePath)
{
    // Check that the file format is matrix market; the only format we can read right now
    // This is not a complete solution, and fails for directories with file names etc...
    // TODO: Should we use boost filesystem?
    std::string strPath( filePath );
    if( strPath.find_last_of( '.' ) != std::string::npos )
    {
        std::string ext = strPath.substr( strPath.find_last_of( '.' ) + 1 );
        if( ext != "mtx" )
        {
            std::cout << "Reading file name error" << std::endl;
            return false;
        }
    }
    else
        return false;

    // Read data from a file on disk into buffers
    // Data is read natively as COO format with the reader
    MatrixMarketReader mm_reader;
    if( mm_reader.MMReadFormat(filePath) )
        return false;

    // JPA: Shouldn't that just be an assertion check? It seems to me that
    // the user have to call clsparseHeaderfromFile before calling this function,
    // otherwise the whole pCsrMatrix will be broken;
    A->nrow = mm_reader.GetNumRows( );
    A->ncol = mm_reader.GetNumCols( );
    A->nnze = mm_reader.GetNumNonZeroes( );

    A->row = (int32_t *)malloc((A->nrow + 1) * sizeof(int32_t));
    A->val = (float *)malloc(A->nnze * sizeof(float));
    A->col = (int32_t *)malloc(A->nnze * sizeof(int32_t));

    if(A->row == NULL || A->col == NULL || A->val == NULL)
    {
        if(A->row == NULL)
            free((void *)A->row);
        if(A->col == NULL)
            free((void *)A->col);
        if(A->val == NULL)
            free((void *)A->val);
        return false;
    }

    //  The following section of code converts the sparse format from COO to CSR
    Coordinate* coords = mm_reader.GetUnsymCoordinates( );

    std::sort( coords, coords + A->nnze, CoordinateCompare );

    int32_t current_row = 1;

    A->row[ 0 ] = 0;

    for(int32_t i = 0; i < A->nnze; i++)
    {
        A->col[ i ] = coords[ i ].y;
        A->val[ i ] = coords[ i ].val;

        while( coords[ i ].x >= current_row )
            A->row[ current_row++ ] = i;
    }

    A->row[ current_row ] = A->nnze;

    while( current_row <= A->nrow )
        A->row[ current_row++ ] = A->nnze;

    return true;
}

uint64_t GetTimeStamp() {
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

void multiply_single(struct sparse_mtx *A, struct dense_mtx *B, struct dense_mtx *C)
{
	// TODO: Implement matrix multiplication with single thread. C=A*B
	C->val = (float *) malloc(sizeof(float) * C->nrow * C->ncol);
	memset(C->val, 0, sizeof(float) * C->nrow * C->ncol);

	for (int r = 0; r < A->nrow; r++) 
		for (int i = A->row[r]; i < A->row[r+1]; i++) 
			for (int j = 0; j < C->ncol; j++)
				C->val[r*C->ncol+j] += A->val[i] * B->val[A->col[i]*B->ncol+j];
}

void multiply_pthread(struct sparse_mtx *A, struct dense_mtx *B, struct dense_mtx *C, int p)
{
	// TODO: Implement matrix multiplication with pthread. C=A*B
	int i, j, k;
	struct args aux[p];
	pthread_t p_threads[p];
	float (*local_sum)[C->ncol];

	C->val = (float *) malloc(sizeof(float) * A->nrow * B->ncol);
	local_sum = (float (*)[C->ncol]) malloc(sizeof(float) * A->nrow * B->ncol);

	memset(C->val, 0, sizeof(float) * C->nrow * C->ncol);
	memset(local_sum, 0, sizeof(float) * C->nrow * C->ncol);

	for (i = 1; i < p; i++) {
		aux[i].A = A;
		aux[i].B = B;
		aux[i].C = C;
		aux[i].p = p;
		aux[i].t_num = i;
		pthread_create(&p_threads[i], NULL, mul_sparse_dense, &aux[i]);
	}

	int t_num = 0;
	int start, num_element;

	start = 0;
	num_element = A->nnze / p;
	if ((t_num + 1) == p)
		num_element += A->nnze % p;
	int row, last_row;
	int cur_val = start;
	row = find_row(A, start);
	last_row = find_row(A, start + num_element - 1);
	for (; row <= last_row; row++) {
		for (i = cur_val; (i < start + num_element) && (i < A->row[row+1]); i++, cur_val++) {
			for (j = 0; j < C->ncol; j++) {
				local_sum[row][j] += A->val[i] * B->val[A->col[i]*B->ncol+j];
			}

		}
	}

	for (i = 1; i < p; i++)
		pthread_join(p_threads[i], NULL);

	for (i = 0; i < C->nrow; i++)
		for (j = 0; j < C->ncol; j++)
			C->val[i*C->ncol+j] += local_sum[i][j];


}



int main(int argc, char **argv)
{
	struct sparse_mtx A;
	int p;
	if(!SCsrMatrixfromFile(&A, argv[1]))
	{
		std::cout << "read failed." << std::endl;
		return 0;
	}

	struct dense_mtx B;
	B.nrow = A.ncol;
	B.ncol = atoi(argv[2]);
	if(B.ncol < 0)
	{
		free(A.row);
		free(A.col);
		free(A.val);
		std::cerr << "Invalid argument for the number of columns of B." << std::endl;
	}
	B.val = (float *)malloc(sizeof(float) * B.nrow * B.ncol);

	srand((unsigned int)time(NULL));
	for(int i = 0; i < B.nrow; i++)
	{
		for(int j = 0; j < B.ncol; j++)
		{
			B.val[B.ncol * i + j] = ((float)rand()/(float)(RAND_MAX)) * ((rand() % 2) ? 1.0f : -1.0f);
		}
	}

	struct dense_mtx C1, C2;
	C1.val = NULL;
	C1.nrow = A.nrow;
	C1.ncol = B.ncol;
	C2.val = NULL;
	C2.nrow = A.nrow;
	C2.ncol = B.ncol;

	p = atoi(argv[3]);

	std::cout << "Single Thread Computation Start" << std::endl;
	uint64_t start = GetTimeStamp();
	multiply_single(&A, &B, &C1);
	uint64_t end = GetTimeStamp();
	std::cout << "Single Thread Computation End: " << end - start  << " us." << std::endl;
	std::cout << "Multi Thread Computation Start" << std::endl;
	start = GetTimeStamp();
	multiply_pthread(&A, &B, &C2, p);
	end = GetTimeStamp();
	std::cout << "Multi Thread Computation End: " << end - start << " us." << std::endl;

	// TODO: Testing Code by comparing C1 and C2
	bool is_correct = true;
	for (int i = 0; i < C1.nrow && is_correct == true; i++)
		for (int j = 0; j < C1.ncol; j++)
			if (!nearly_equal(C1.val[i*C1.ncol+j], C2.val[i*C2.ncol+j], 0)) {
				printf("%lf %lf\n", C1.val[i*C1.ncol+j], C2.val[i*C2.ncol+j]);
				is_correct = false;
				break;
			}

	printf("%s\n", is_correct ? "correct" : "wrong");

	free(A.row);
	free(A.col);
	free(A.val);
	free(B.val);
	if(C1.val != NULL)
		free(C1.val);
	if(C2.val != NULL)
		free(C2.val);

	return 0;
}

void *mul_sparse_dense(void *arg)
{
	int i, j, k;
	struct args *aux = (struct args *) arg;
	struct sparse_mtx *A = aux->A;
	struct dense_mtx *B = aux->B;
	struct dense_mtx *C = aux->C;
	int p = aux->p;
	int t_num = aux->t_num;
	float (*local_sum)[C->ncol];

	local_sum = (float (*)[C->ncol]) malloc(sizeof(float) * A->nrow * B->ncol);
	memset(local_sum, 0, sizeof(float) * C->nrow * C->ncol);

	int start, num_element;

	num_element = A->nnze / p;
	start = num_element * t_num;
	if ((t_num + 1) == p)
		num_element += A->nnze % p;
	int row, last_row;
	int cur_val = start;
	row = find_row(A, start);
	last_row = find_row(A, start + num_element - 1);
	for (; row <= last_row; row++) {
		for (i = cur_val; (i < start + num_element) && (i < A->row[row+1]); i++, cur_val++) {
			for (j = 0; j < C->ncol; j++) {
				local_sum[row][j] += A->val[i] * B->val[A->col[i]*B->ncol+j];
			}

		}
	}

	pthread_mutex_lock(&C_lock);
	for (i = 0; i < C->nrow; i++)
		for (j = 0; j < C->ncol; j++)
			C->val[i*C->ncol+j] += local_sum[i][j];
	pthread_mutex_unlock(&C_lock);

	return NULL;
}
int find_row(struct sparse_mtx *A, int val_i)
{
	int i;

	for (i = 0; i < A->nrow; i++)
		if ((A->row[i] <= val_i) && (val_i < A->row[i+1]))
			return i;
	return -1;
}

bool nearly_equal(float a, float b, float epsilon) 
{
	float abs_a = abs(a);
	float abs_b = abs(b);
	float diff = abs(a - b);
	float float_min_normal;
	float float_max;
	long temp = 1;
	float_min_normal = *((float *) &temp);
	float_max = numeric_limits<float>::max();

	if (epsilon == 0)
		epsilon = 0.00001;

	if (a == b) { // shortcut, handles infinities
		return true;
	} else if (a == 0 || b == 0 || diff < float_min_normal) {
		// a or b is zero or both are extremely close to it
		// relative error is less meaningful here
		return diff < (epsilon * float_min_normal);
	} else { // use relative error
		return diff / min((abs_a + abs_b), float_max) < epsilon;
	}
}
