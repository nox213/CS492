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

void swap_row(int n, double *arr, double *b, int x, int y);
void single_elimination(int n, double *a, double *b);
void back_substitution(int n, double *a, double *b, double *x);
int maxloc(int start, int n, double *a);
bool nearly_equal(double a, double b, double epsilon); 

uint64_t GetTimeStamp() {
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

int main(int argc, char *argv[])
{
	int n, p;
	int i, j, k;
	uint64_t begin, end;
	uint64_t elapsed_s, elapsed_p;

	if (argc < 3) {
		fprintf(stderr, "gaussian.out n p\n");
		return 0;
	}

	n = strtol(argv[1], NULL, 10);
	p = strtol(argv[2], NULL, 10);

	double (*a)[n], (*b), (*x), (*serial_a)[n], (*serial_b), (*serial_x), (*A)[n], (*B);

	if ((a = (double (*)[n]) malloc(sizeof(double) * (n * n))) == NULL) {
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

	memset(a, 0, sizeof(double) * (n * n));
	memset(b, 0, sizeof(double) * n);
	memset(x, 0, sizeof(double) * n);
	memset(serial_a, 0, sizeof(double) * (n * n));
	memset(serial_b, 0, sizeof(double) * n);
	memset(serial_x, 0, sizeof(double) * n);

	srand48(time(NULL));
	for (i = 0; i < n; i++) {
		B[i] = serial_b[i] = b[i] = drand48();
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

	for (j = 0; j < n - 1; j++) {
		int row_max = maxloc(j, n, (double *) a);
		swap_row(n, (double *) a, (double *) b, j, row_max); 
		for (i = j + 1; i < n; i++) {
			double m = a[i][j] / a[j][j];
			for (k = j; k < n; k++) 
				a[i][k] -= m * a[j][k];
			b[i] -= m * b[j];
		}
	}
	back_substitution(n, (double *) a, (double *) b, x);
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

int maxloc(int start, int n, double *a)
{
	int i, row_max, col = start;
	double cur_max;

	row_max = start;
	cur_max = fabs(a[start * n + col]);
	for (i = start + 1; i < n; i++)
		if (cur_max < fabs(a[i * n + col])) {
			cur_max = fabs(a[i * n + col]);
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
