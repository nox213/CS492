#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include <stdint.h>

#define NANO_TO_SEC(x) ((double) ((x) / 10e9))
#define NANO_TO_MILI(x) ((long) ((x) / 10e6))
#define SEC_TO_NANO(x) (((x) * 1000000000))
#define SEC_TO_MILI(x) ((long) ((x) * 10e3))

void print_array(int n, double arr[][n]);
void swap_row(int n, double arr[][n], double b[n], int x, int y);
void single_elimination(int n, double a[][n], double b[n]);
void back_substitution(int n, double a[][n], double b[n], double x[n]);
int maxloc(int start, int n, double (*a)[n]);
bool nearly_equal(double a, double b, double epsilon); 

uint64_t GetTimeStamp() {
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

static inline double min(const double a, const double b)
{
	return a < b ? a : b;
}

static inline double max(const double a, const double b)
{
	return a > b ? a : b;
}

int main(int argc, char *argv[])
{
	int n, p;
	int i, j, k;
	uint64_t begin, end, begin_b, end_b;
	long elapsed_s, elapsed_p, elapsed_b;

	if (argc < 3) {
		fprintf(stderr, "gaussian.out n p\n");
		return 0;
	}

	n = strtol(argv[1], NULL, 10);
	p = strtol(argv[2], NULL, 10);

	double (*a)[n], (*b), (*x), (*serial_a)[n], (*serial_b), (*serial_x), (*A)[n], (*B);

	if ((a = malloc(sizeof(double) * (n * n))) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}
	if ((b = malloc(sizeof(double) * n)) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}
	if ((x = malloc(sizeof(double) * n)) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}

	if ((serial_a = malloc(sizeof(double) * (n * n))) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}

	if ((serial_b = malloc(sizeof(double) * n)) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}
	if ((serial_x = malloc(sizeof(double) * n)) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}
	if ((A = malloc(sizeof(double) * (n * n))) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}
	if ((B = malloc(sizeof(double) * n)) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}

	memset(a, 0, sizeof(double) * (n * n));
	memset(b, 0, sizeof(double) * n);
	memset(x, 0, sizeof(double) * n);
	memset(serial_a, 0, sizeof(double) * (n * n));
	memset(serial_b, 0, sizeof(double) * n);
	memset(serial_x, 0, sizeof(double) * n);

	for (i = 0; i < n; i++) {
		B[i] = serial_b[i] = b[i] = drand48();
		for (j = 0; j < n; j++) {
			A[i][j] = serial_a[i][j] = a[i][j] = drand48();
		}
	}

	printf("Single thread computaion start\n");
	begin = GetTimeStamp();
	single_elimination(n, serial_a, serial_b);
	back_substitution(n, serial_a, serial_b, serial_x);
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

	omp_set_num_threads(p);
	for (j = 0; j < n - 1; j++) {
		int row_max = maxloc(j, n, a);
		swap_row(n, a, b, j, row_max); 
		#pragma omp parallel for private(i, k) schedule(static)
		for (i = j + 1; i < n; i++) {
			double m = a[i][j] / a[j][j];
			for (k = j; k < n; k++) 
				a[i][k] -= m * a[j][k];
			b[i] -= m * b[j];
		}
	}
	begin_b = GetTimeStamp();
	back_substitution(n, a, b, x);
	end_b = GetTimeStamp();

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

void print_array(int n, double arr[][n])
{
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) 
			printf("%g ", arr[i][j]);
		printf("\n");
	}
}



bool nearly_equal(double a, double b, double epsilon) 
{
	double abs_a = fabs(a);
	double abs_b = fabs(b);
	double diff = fabs(a - b);
	double double_min_normal;
	long temp = 1;

	double_min_normal = *((double *) &temp);

	if (epsilon == 0)
		epsilon = 0.00001;

	if (a == b) { // shortcut, handles infinities
		return true;
	} else if (a == 0 || b == 0 || diff < double_min_normal) {
		// a or b is zero or both are extremely close to it
		// relative error is less meaningful here
		return diff < (epsilon * double_min_normal);
	} else { // use relative error
		return diff / min((abs_a + abs_b), DBL_MAX) < epsilon;
	}
}

int maxloc(int start, int n, double (*a)[n])
{
	int i, row_max, col = start;
	double cur_max;

	row_max = start;
	cur_max = fabs(a[start][col]);
	for (i = start + 1; i < n; i++)
		if (cur_max < fabs(a[i][col])) {
			cur_max = fabs(a[i][col]);
			row_max = i;
		}

	return row_max;
}

void swap_row(int n, double arr[][n], double b[n], int x, int y)
{
	double *temp, t;

	if ((temp = malloc(sizeof(double) * n)) == NULL) {
		fprintf(stderr, "malloc error in %d\n", __LINE__);
		exit(EXIT_FAILURE);
	}

	memcpy(temp, arr[x], sizeof(arr[x]));
	memcpy(arr[x], arr[y], sizeof(arr[x]));
	memcpy(arr[y], temp, sizeof(arr[x]));

	t = b[x] ;
	b[x] = b[y];
	b[y] = t;
}

void back_substitution(int n, double a[][n], double b[n], double x[n])
{
	int i, j;

	for (i = n - 1; i >= 0; i--) {
		x[i] = b[i];
		for (j = i + 1; j < n; j++) 
			x[i] -= a[i][j] * x[j];
		x[i] /= a[i][i];
	}
}

void single_elimination(int n, double a[][n], double b[n])
{
	int i, j, k;

	for (j = 0; j < n - 1; j++) {
		int row_max = maxloc(j, n, a);
		swap_row(n, a, b, j, row_max); 
		for (i = j + 1; i < n; i++) {
			double m = a[i][j] / a[j][j];
			for (k = j; k < n; k++) 
				a[i][k] -= m * a[j][k];
			b[i] -= m * b[j];
		}
	}
}
