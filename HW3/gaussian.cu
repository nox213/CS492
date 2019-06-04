#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <stdint.h>

#define SEC_TO_NANO(x) (((x) * 1000000000L))
#define TILE_WIDTH 16
#define BLOCK_SIZE 32

void swap_row(int n, double *arr, double *b, int x, int y);
void back_substitution(int n, double *a, double *b, double *x);
int maxloc(int start, int n, double *a);
bool nearly_equal(double a, double b, double epsilon); 
__global__ void copy_col_kernel(double *temp, double *a, int n, int j);
__global__ void swap_row_1(double *a, double *temp, int n, int row1);
__global__ void swap_row_2(double *a, double *temp, int n, int row1, int row2);
__global__ void swap_row_3(double *a, double *temp, int n, int row2);
__global__ void make_mul(double *a, double *m, int n, int j);
__global__ void gaussian_elimination(double *a, double *m, int n, int j);
__global__ void pivot_to_one(double *a, int n, double *temp, int pivot);

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
		fprintf(stderr, "%d", argc);
		fprintf(stderr, "gaussian.out n\n");
		return 0;
	}
	n = strtol(argv[1], NULL, 10);

	double (*a)[n+1], (*b), (*x), (*serial_a)[n], (*serial_b), (*serial_x), (*A)[n], (*B);

	if ((a = (double (*)[n+1]) malloc(sizeof(double) * ((n + 1) * n))) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}
	if ((b = (double *) malloc(sizeof(double) * n)) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}
	if ((x = (double *) malloc(sizeof(double) * n)) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}

	if ((serial_a = (double (*)[n]) malloc(sizeof(double) * (n * n))) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}

	if ((serial_b = (double *) malloc(sizeof(double) * n)) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}
	if ((serial_x = (double *) malloc(sizeof(double) * n)) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}
	if ((A = (double (*)[n]) malloc(sizeof(double) * (n * n))) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}
	if ((B = (double *) malloc(sizeof(double) * n)) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}

	memset(a, 0, sizeof(double) * ((n + 1) * n));
	memset(b, 0, sizeof(double) * n);
	memset(x, 0, sizeof(double) * n);
	memset(serial_a, 0, sizeof(double) * (n * n));
	memset(serial_b, 0, sizeof(double) * n);
	memset(serial_x, 0, sizeof(double) * n);

	srand48(time(NULL));
	for (i = 0; i < n; i++) {
		a[i][n] = B[i] = serial_b[i] = b[i] = drand48();
		for (j = 0; j < n; j++) {
			A[i][j] = serial_a[i][j] = a[i][j] = drand48();
		}
	}

	printf("Single thread computaion start\n");
	begin = GetTimeStamp();
	for (j = 0; j < n - 1; j++) {
		int row_max = maxloc(j, n, (double *) serial_a);
		swap_row(n, (double *) serial_a, (double *) serial_b, j, row_max); 
		for (i = j + 1; i < n; i++) {
			double m = serial_a[i][j] / serial_a[j][j];
			for (k = j; k < n; k++) 
				serial_a[i][k] -= m * serial_a[j][k];
			serial_b[i] -= m * serial_b[j];
		}
	}
	back_substitution(n, (double *) serial_a, (double *) serial_b, (double *) serial_x);
	end = GetTimeStamp();
	elapsed_s = end - begin;
	printf("Single thread computaion end\n");
	printf("elapsed time: %ld (usec)\n", elapsed_s);

	double l2norm = 0;
	for (i = 0; i < n; i++) {
		double ax = 0;
		for (j = 0; j < n; j++)
			ax += A[i][j] * serial_x[j];
		l2norm += (ax - B[i]) * (ax - B[i]);

	}
	printf("l2norm of single: %g\n", sqrt(l2norm));


	printf("Multi thread computaion start\n");
	begin = GetTimeStamp();

	double *device_a, *device_b, *device_temp, *col, *device_m;
	dim3 dimGrid(n / TILE_WIDTH + 1, (n + 1) / TILE_WIDTH + 1, 1);
	dim3 dimBlock(TILE_WIDTH, TILE_WIDTH, 1);

	if (cudaErrorMemoryAllocation == cudaMalloc(&device_a, sizeof(double) * n * (n + 1)))
		fprintf(stderr, "error\n");
	if (cudaErrorMemoryAllocation == cudaMalloc(&device_b, sizeof(double) * n))
		fprintf(stderr, "error\n");
	if (cudaErrorMemoryAllocation == cudaMalloc(&device_temp, sizeof(double) * (n + 1)))
		fprintf(stderr, "error\n");
	if (cudaErrorMemoryAllocation == cudaMalloc(&device_m, sizeof(double) * n))
		fprintf(stderr, "error\n");

	col = (double *) malloc(sizeof(double) * n);
	cudaMemcpy(device_a, (void **) a, sizeof(double) * (n + 1) * n, cudaMemcpyHostToDevice);
	for (j = 0; j < n; j++) {
		copy_col_kernel<<<n / BLOCK_SIZE + 1, BLOCK_SIZE>>>(device_temp, device_a, n, j);
		cudaMemcpy(col, device_temp, sizeof(double) * n, cudaMemcpyDeviceToHost);
		int row_max = maxloc(j, n, col);
		swap_row_1<<<(n + 1) / BLOCK_SIZE + 1, BLOCK_SIZE>>>(device_a, device_temp, n, row_max);
		swap_row_2<<<(n + 1) / BLOCK_SIZE + 1, BLOCK_SIZE>>>(device_a, device_temp, n, row_max, j);
		swap_row_3<<<(n + 1) / BLOCK_SIZE + 1, BLOCK_SIZE>>>(device_a, device_temp, n, j);
		copy_col_kernel<<<n / BLOCK_SIZE + 1, BLOCK_SIZE>>>(device_temp, device_a, n, j);
		pivot_to_one<<<(n + 1) / BLOCK_SIZE + 1, BLOCK_SIZE>>>(device_a, n, device_temp, j);
		make_mul<<<n / BLOCK_SIZE + 1, BLOCK_SIZE>>>(device_a, device_m, n, j);
		gaussian_elimination<<<dimGrid, dimBlock>>>(device_a, device_m, n, j);
	}
	copy_col_kernel<<<n / BLOCK_SIZE + 1, BLOCK_SIZE>>>(device_b, device_a, n, n);
	cudaMemcpy(x, device_b, sizeof(double) * n, cudaMemcpyDeviceToHost);
	end = GetTimeStamp();

	elapsed_p = end - begin;
	printf("Multi thread computaion end\n");
	printf("elapsed time: %ld (usec)\n", elapsed_p);

	l2norm = 0;
	for (i = 0; i < n; i++) {
		double ax = 0;
		for (j = 0; j < n; j++)
			ax += A[i][j] * x[j];
		l2norm += (ax - B[i]) * (ax - B[i]);

	}
	printf("l2norm of multi: %g\n", sqrt(l2norm));

	printf("Speed up is %g\n", (double) elapsed_s / elapsed_p);

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

int maxloc(int start, int n, double *col)
{
	int i, row_max;
	double cur_max;

	row_max = start;
	cur_max = fabs(col[start]);
	for (i = start + 1; i < n; i++)
		if (cur_max < fabs(col[i])) {
			cur_max = fabs(col[i]);
			row_max = i;
		}

	return row_max;
}

void swap_row(int n, double *arr, double *b, int x, int y)
{
	double *temp, t;
	int size = sizeof(double) * n;

	if ((temp = (double *) malloc(sizeof(double) * n)) == NULL) {
		fprintf(stderr, "malloc error in %d\n", __LINE__);
		exit(EXIT_FAILURE);
	}

	memcpy(temp, arr + x * n, size);
	memcpy(arr + x * n, arr + y * n, size);
	memcpy(arr + y * n, temp, size);

	t = b[x] ;
	b[x] = b[y];
	b[y] = t;
}

void back_substitution(int n, double *a, double *b, double *x)
{
	int i, j;

	for (i = n - 1; i >= 0; i--) {
		x[i] = b[i];
		for (j = i + 1; j < n; j++) 
			x[i] -= a[i * n + j] * x[j];
		x[i] /= a[i * n + i];
	}
}

__global__ void copy_col_kernel(double *temp, double *a, int n, int j)
{
	int bx = blockIdx.x;
	int tx = threadIdx.x;

	int row = bx * blockDim.x + tx;

	if (row < n)
		temp[row] = a[row * (n + 1) + j];
}

__global__ void swap_row_1(double *a, double *temp, int n, int row1)
{
	int bx = blockIdx.x;
	int tx = threadIdx.x;

	int col  = bx * blockDim.x + tx;

	if (col < n + 1)
		temp[col] = a[row1 * (n + 1) + col];
}

__global__ void swap_row_2(double *a, double *temp, int n, int row1, int row2)
{
	int bx = blockIdx.x;
	int tx = threadIdx.x;

	int col  = bx * blockDim.x + tx;

	if (col < n + 1)
		a[row1 * (n + 1) + col] = a[row2 * (n + 1) + col];
}

__global__ void swap_row_3(double *a, double *temp, int n, int row2)
{
	int bx = blockIdx.x;
	int tx = threadIdx.x;

	int col  = bx * blockDim.x + tx;

	if (col < n + 1)
		a[row2 * (n + 1) + col] = temp[col];
}

__global__ void make_mul(double *a, double *m, int n, int j)
{
	int bx = blockIdx.x;
	int tx = threadIdx.x;

	int row  = bx * blockDim.x + tx;

	if (row < n)
		m[row] = a[row * (n + 1) + j] / a[j * (n + 1) + j];
}

__global__ void gaussian_elimination(double *a, double *m, int n, int j)
{
	//	__shared__ double sub_m[TILE_WIDTH];
	//	__shared__ double 

	int bx = blockIdx.x;
	int by = blockIdx.y;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int row = by * TILE_WIDTH + ty; 
	int col = bx * TILE_WIDTH + tx;

	if (row != j && row < n && col < (n + 1) && col >= j)
		a[row * (n + 1) + col] -= m[row] * a[j * (n + 1) + col];
}

__global__ void pivot_to_one(double *a, int n, double* temp, int pivot)
{
	int bx = blockIdx.x;
	int tx = threadIdx.x;

	int col  = bx * blockDim.x + tx;

	if (col < n + 1)
		a[pivot * (n + 1) + col] /= temp[pivot];
}
