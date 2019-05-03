#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#define NANO_TO_SEC(x) ((double) ((x) / 10e9))
#define NANO_TO_MILI(x) ((long) ((x) / 10e6))
#define SEC_TO_NANO(x) ((long) ((x) * 10e9))
#define SEC_TO_MILI(x) ((long) ((x) * 10e3))

void mul_matrix_single(int n, double (*a)[n], double (*b)[n], double (*answer)[n]);
bool nearly_equal(double a, double b, double epsilon); 
void transepose(int n, double (*m)[n]);

static inline double min(const double a, const double b)
{
	return a < b ? a : b;
}

int main(int argc, char *argv[])
{
	int n, p;
	int i, j, k;
	struct timespec begin, end;
	long elapsed_s, elapsed_p;

	if (argc < 3) {
		fprintf(stderr, "dense n p\n");
		return 0;
	}

	n = strtol(argv[1], NULL, 10);
	p = strtol(argv[2], NULL, 10);

	double (*a)[n], (*b)[n], (*c)[n], (*answer)[n];

	if ((a = malloc(sizeof(double) * (n * n))) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}
	if ((b = malloc(sizeof(double) * (n * n))) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}
	if ((c = malloc(sizeof(double) * (n * n))) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}

	if ((answer = malloc(sizeof(double) * (n * n))) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}

	memset(a, 0, sizeof(double) * (n * n));
	memset(b, 0, sizeof(double) * (n * n));
	memset(c, 0, sizeof(double) * (n * n));
	memset(answer, 0, sizeof(double) * (n * n));

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			a[i][j] = drand48();
			b[i][j] = drand48();
		}
	}

	printf("Single thread computaion start\n");
	clock_gettime(CLOCK_REALTIME, &begin);
	mul_matrix_single(n, a, b, answer);
	clock_gettime(CLOCK_REALTIME, &end);
	elapsed_s = SEC_TO_NANO(end.tv_sec - begin.tv_sec) + ((end.tv_nsec - begin.tv_nsec));
	printf("Single thread computaion end\n");
	printf("elapsed time: %lf (sec)\n", NANO_TO_SEC(elapsed_s));

	double sum;
	printf("Multi thread computaion start\n");
	clock_gettime(CLOCK_REALTIME, &begin);
	transepose(n, b);
	omp_set_num_threads(p);
#pragma omp parallel for shared(a, b, c, n) private(i, j, k, sum) schedule(static)
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				sum = 0;
				for (k = 0; k < n; k++) {
					sum += a[i][k] * b[j][k];
				}
				c[i][j] = sum;
			}
		}
	clock_gettime(CLOCK_REALTIME, &end);

	elapsed_p = SEC_TO_NANO(end.tv_sec - begin.tv_sec) + ((end.tv_nsec - begin.tv_nsec));
	printf("Multi thread computaion end\n");
	printf("elapsed time: %lf (sec)\n", NANO_TO_SEC(elapsed_p));

	printf("\nSpeed up is %g\n", (double) elapsed_s / elapsed_p);

	bool is_correct = true;
	for (i = 0; i < n && is_correct == true; i++)
		for (j = 0; j < n; j++)
			if (!nearly_equal(c[i][j], answer[i][j], 0)) {
				printf("%lf %lf\n", c[i][j], answer[i][j]);
				is_correct = false;
				break;
			}

	printf("%s\n", is_correct ? "correct" : "wrong");

	return 0;
}


void mul_matrix_single(int n, double (*a)[n], double (*b)[n], double (*answer)[n])
{
	int i, j, k;

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			for (k = 0; k < n; k++)
				answer[i][j] += a[i][k] * b[k][j];
}

void transepose(int n, double (*m)[n])
{
	int i, j;
	double (*temp)[n];

	temp = malloc(sizeof(double) * n * n);

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			temp[j][i] = m[i][j];

	memcpy(m, temp, sizeof(double) * n * n);
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
	} else if (a == 0 || b == 0 || diff < DBL_MIN) {
		// a or b is zero or both are extremely close to it
		// relative error is less meaningful here
		return diff < (epsilon * DBL_MIN);
	} else { // use relative error
		return diff / min((abs_a + abs_b), DBL_MAX) < epsilon;
	}
}
