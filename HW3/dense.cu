#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <stdint.h>
#include <time.h>

#define TILE_WIDTH 32

bool nearly_equal(double a, double b, double epsilon); 
__global__ void MatrixMulKernel(double *M, double *N, double *P, int Width);

uint64_t GetTimeStamp() {
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

int main(int argc, char *argv[])
{
	int n;
	int i, j, k;
	uint64_t begin, end;
	uint64_t elapsed_s, elapsed_p;

	if (argc < 2) {
		fprintf(stderr, "dense n\n");
		return 0;
	}

	n = strtol(argv[1], NULL, 10);

	double (*a)[n], (*b)[n], (*c)[n], (*answer)[n];
	int size = sizeof(double) * n * n;

	if ((a = (double (*)[n]) malloc(size)) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}
	if ((b = (double (*)[n]) malloc(size)) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}
	if ((c = (double (*)[n]) malloc(size)) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}

	if ((answer = (double (*)[n]) malloc(size)) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}

	memset(a, 0, size);
	memset(b, 0, size);
	memset(c, 0, size);
	memset(answer, 0, size);

	srand48(time(NULL));
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			a[i][j] = drand48();
			b[i][j] = drand48();
		}
	}

	int grid_size = n / TILE_WIDTH;
	if (n % TILE_WIDTH)
		grid_size++;

	double *A, *B, *C;
	dim3 dimGrid(grid_size, grid_size, 1);
	dim3 dimBlock(TILE_WIDTH, TILE_WIDTH, 1);

	if (cudaErrorMemoryAllocation == cudaMalloc((void **) &A, size))
		fprintf(stderr, "error\n");
	if (cudaErrorMemoryAllocation == cudaMalloc((void **) &B, size))
		fprintf(stderr, "error\n");
	if (cudaErrorMemoryAllocation == cudaMalloc((void **) &C, size))
		fprintf(stderr, "error\n");

	printf("Single thread computaion start\n");
	begin = GetTimeStamp();
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			for (k = 0; k < n; k++)
				answer[i][j] += a[i][k] * b[k][j];
	end = GetTimeStamp();
	elapsed_s = end - begin;
	printf("Single thread computaion end\n");
	printf("elapsed time: %ld (usec)\n", elapsed_s);

	printf("Multi thread computaion start\n");
	cudaMemcpy(A, (void **) a, size, cudaMemcpyHostToDevice);
	cudaMemcpy(B, (void **) b, size, cudaMemcpyHostToDevice);

	begin = GetTimeStamp();
	MatrixMulKernel<<<dimGrid, dimBlock>>>(A, B, C, n);
	cudaMemcpy(c, C, size, cudaMemcpyDeviceToHost);
	end = GetTimeStamp();

	elapsed_p = end - begin;
	printf("Multi thread computaion end\n");
	printf("elapsed time: %ld (usec)\n", elapsed_p);

	printf("\nSpeed up is %g\n", (double) elapsed_s / elapsed_p);

	bool is_correct = true;
	for (i = 0; i < n && is_correct == true; i++)
		for (j = 0; j < n; j++)
			if (!nearly_equal(c[i][j], answer[i][j], 0)) {
				printf("i: %d j : %d\n", i, j);
				printf("%g %g\n", c[i][j], answer[i][j]);
				is_correct = false;
				break;
			}

	printf("%s\n", is_correct ? "correct" : "wrong");
	cudaFree(A);
	cudaFree(B);
	cudaFree(C);

	return 0;
}

bool nearly_equal(double a, double b, double epsilon) 
{
	double abs_a = fabs(a);
	double abs_b = fabs(b);
	double diff = fabs(a - b);

	if (epsilon == 0)
		epsilon = 0.00001;

	if (a == b) { // shortcut, handles infinities
		return true;
	} else if (a == 0 || b == 0 || diff < DBL_MIN) {
		// a or b is zero or both are extremely close to it
		// relative error is less meaningful here
		return diff < (epsilon * DBL_MIN);
	} else { // use relative error
		return diff / min((abs_a + abs_b), DBL_MAX) < epsilon;
	}
}

__global__ void MatrixMulKernel(double *M, double *N, double *P, int Width)
{
	__shared__ double subTileM[TILE_WIDTH][TILE_WIDTH];
	__shared__ double subTileN[TILE_WIDTH][TILE_WIDTH];

	int bx = blockIdx.x;
	int by = blockIdx.y;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int Row = by * TILE_WIDTH + ty;
	int Col = bx * TILE_WIDTH + tx;
	double Pvalue = 0;
	int num_tile = Width / TILE_WIDTH;
	if (Width % TILE_WIDTH)
		num_tile++;

	for (int m = 0; m < num_tile; m++) {
		if (m * TILE_WIDTH + tx < Width ) 
			subTileM[ty][tx] = M[Row*Width+m*TILE_WIDTH+tx];
		else
			subTileM[ty][tx] = 0;
		if (m * TILE_WIDTH + ty < Width)
			subTileN[ty][tx] = N[(m*TILE_WIDTH+ty)*Width+Col];
		else 
			subTileN[ty][tx] = 0;
		__syncthreads();
		for (int k = 0; k < TILE_WIDTH; k++)
			Pvalue += subTileM[ty][k] * subTileN[k][tx];
		__syncthreads();
	}
	if (Row < Width && Col < Width)
		P[Row*Width+Col] = Pvalue;
}

