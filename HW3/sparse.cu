#include "mmreader.hpp"
#include <time.h>
#include <iostream>
#include <sys/time.h>
#include <unistd.h>
#include <algorithm>
#include <cfloat>
#include <omp.h>

#define TILE_WIDTH 16

using namespace std;

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

	for (int32_t i = 0; i < A->nnze; i++)
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

void multiply_single(struct sparse_mtx *A, struct dense_mtx *B, struct dense_mtx *C)
{
	C->nrow = A->nrow;
	C->ncol = B->ncol;
	C->val = (float *)malloc(C->nrow * C->ncol * sizeof(float));

	if(C->val == NULL)
		return;

	memset(C->val, 0, sizeof(float) * C->nrow * C->ncol);
	// TODO: Implement matrix multiplication with single thread. C=A*B
	for(int32_t i = 0; i < A->nrow; i++)
	{
		int32_t A_col_start = A->row[i];
		int32_t A_col_stop = A->row[i + 1];

		for(int32_t j = A_col_start; j < A_col_stop; j++)
		{
			int32_t B_row = A->col[j];

			for(int32_t k = 0; k < B->ncol; k++) 
				C->val[i * C->ncol + k] += A->val[j] * B->val[B_row * B->ncol + k];
		}
	}
}

__global__ void multiply_kernel(int32_t *a_row, int32_t *a_col, float *a_val, float *b_val, float *c_val, 
		int a_nrow, int a_ncol, int b_nrow, int b_ncol)
{
	int bx = blockIdx.x;
	int by = blockIdx.y;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int row = by * blockDim.y + ty;
	int col = bx * blockDim.x + tx;
	float temp = 0;

	if (row < a_nrow && col < b_ncol)
		for (int i = a_row[row]; i < a_row[row+1]; i++) 
			temp += a_val[i] * b_val[a_col[i]*b_ncol+col];
	c_val[row*b_ncol+col] = temp;
}

uint64_t GetTimeStamp() {
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

int main(int argc, char **argv)
{
	struct sparse_mtx A;

	if(!SCsrMatrixfromFile(&A, argv[1]))
	{
		std::cout << "read failed." << std::endl;
		return 0;
	}
	std::cout << "Matrix: " << argv[1] << std::endl;

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

	uint64_t time_s, time_p;

	std::cout << "Single Thread Computation Start" << std::endl;
	uint64_t start = GetTimeStamp();
	multiply_single(&A, &B, &C1);
	uint64_t end = GetTimeStamp();
	std::cout << "Single Thread Computation End: " << end - start  << " us." << std::endl;
	time_s = end - start;

	int32_t *a_row;
	int32_t *a_col;
	float *a_val;
	float *b_val;
	float *c_val;
	int grid_x, grid_y;

	grid_x = C2.nrow / TILE_WIDTH;
	if (C2.nrow % TILE_WIDTH)
		grid_x++;
	grid_y = C2.ncol / TILE_WIDTH;
	if (C2.ncol % TILE_WIDTH)
		grid_y++;
	dim3 dimGrid(grid_x, grid_y, 1);
	dim3 dimBlock(TILE_WIDTH, TILE_WIDTH, 1);

	cudaMalloc(&a_row, sizeof(int32_t) * A.nrow);
	cudaMalloc(&a_col, sizeof(int32_t) * A.ncol);
	cudaMalloc(&a_val, sizeof(float) * A.nnze);
	cudaMalloc(&b_val, sizeof(float) * B.nrow * B.ncol);
	cudaMalloc(&c_val, sizeof(float) * C2.nrow * C2.ncol);
	C2.val = (float *) malloc(sizeof(float) * C2.nrow * C2.ncol);
	memset(C2.val, 0, sizeof(float)  * C2.nrow * C2.ncol);

	std::cout << "Cuda Computation Start" << std::endl;
	start = GetTimeStamp();

	cudaMemcpy(a_row, A.row, sizeof(int32_t) * A.nrow, cudaMemcpyHostToDevice);
	cudaMemcpy(a_col, A.col, sizeof(int32_t) * A.ncol, cudaMemcpyHostToDevice);
	cudaMemcpy(a_val, A.val, sizeof(float) * A.nnze, cudaMemcpyHostToDevice);
	cudaMemcpy(b_val, B.val, sizeof(float) * B.nrow * B.ncol, cudaMemcpyHostToDevice);

	multiply_kernel<<<dimGrid, dimBlock>>>(a_row, a_col, a_val, b_val, c_val, A.nrow, A.ncol,
			B.nrow, B.ncol);
	cudaMemcpy(C2.val, c_val, sizeof(float) * C2.nrow * C2.ncol, cudaMemcpyDeviceToHost);

	end = GetTimeStamp();
	std::cout << "Cuda Computation End: " << end - start << " us." << std::endl << std::endl;
	time_p = end - start;

	printf("Speed up is %g\n", (double) time_s / time_p);

	// TODO: Testing Code by comparing C1 and C2
	bool is_correct = true;
	for (int i = 0; i < C1.nrow && is_correct == true; i++)
		for (int j = 0; j < C1.ncol; j++)
			if (!nearly_equal(C1.val[i*C1.ncol+j], C2.val[i*C2.ncol+j], 0)) {
				printf("%g %g\n", C1.val[i*C1.ncol+j], C2.val[i*C2.ncol+j]);
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

	cudaFree(a_row);
	cudaFree(a_col);
	cudaFree(a_val);
	cudaFree(b_val);
	cudaFree(c_val);

	return 0;
}

bool nearly_equal(float a, float b, float epsilon) 
{
	float abs_a = abs(a);
	float abs_b = abs(b);
	float diff = abs(a - b);
	float float_min_normal;
	long temp = 1;

	float_min_normal = *((float *) &temp);

	if (epsilon == 0)
		epsilon = 0.00001f;

	if (a == b) { // shortcut, handles infinities
		return true;
	} else if (a == 0 || b == 0 || diff < FLT_MIN) {
		// a or b is zero or both are extremely close to it
		// relative error is less meaningful here
		return diff < (epsilon * float_min_normal);
	} else { // use relative error
		return diff / min((abs_a + abs_b), FLT_MAX) < epsilon;
	}
}
