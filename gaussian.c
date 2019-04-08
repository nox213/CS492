#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>

#define NANO_TO_SEC(x) ((double) ((x) / 10e9))
#define NANO_TO_MILI(x) ((long) ((x) / 10e6))
#define SEC_TO_NANO(x) ((long) ((x) * 10e9))
#define SEC_TO_MILI(x) ((long) ((x) * 10e3))

struct args {
	double **a, *b;
	int n;
	int p;
	int t_num;
};

bool is_ready;

pthread_barrier_t barrier_thread;

void print_array(int n, double arr[][n]);
void swap_row(int n, double arr[][n], double b[n], int x, int y);
void *elimination(void *arg);
void back_substitution(int n, double a[][n], double b[n], double x[n]);
int maxloc(int start, int n, double (*a)[n]);
bool nearly_equal(double a, double b, double epsilon); 

static inline double (const double a, const double b)
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
	struct timespec begin, end;
	long elapsed;

	if (argc < 3) {
		fprintf(stderr, "gaussian.out n p\n");
		return 0;
	}

	n = strtol(argv[1], NULL, 10);
	p = strtol(argv[2], NULL, 10);

	double (*a)[n], (*b), (*x), (*serial)[n], (*serial_b), (*serial_x);
	pthread_t p_threads[p];
	struct args aux[p];

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

	if ((serial = malloc(sizeof(double) * (n * n))) == NULL) {
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

	memset(p_threads, 0, sizeof(p_threads));
	memset(a, 0, sizeof(double) * (n * n));
	memset(b, 0, sizeof(double) * n);
	memset(x, 0, sizeof(double) * n);
	memset(serial, 0, sizeof(double) * (n * n));
	memset(serial_b, 0, sizeof(double) * n);
	memset(serial_x, 0, sizeof(double) * n);

	for (i = 0; i < n; i++) {
		serial_b[i] = b[i] = drand48();
		for (j = 0; j < n; j++) {
			serial[i][j] = a[i][j] = drand48();
		}
	}

	pthread_barrier_init(&barrier_thread, NULL, p);

	clock_gettime(CLOCK_REALTIME, &begin);

	is_ready = false;
	for (i = 1; i < p; i++) {
		aux[i].a = (double **) a;
		aux[i].b =  b;
		aux[i].n = n;
		aux[i].p = p;
		aux[i].t_num = i + 1;
		pthread_create(&p_threads[i], NULL, elimination, (void *) &aux[i]);
	}

	int t_num = 1;
	int start, block_size;
	for (j = 0; j < n - 1; j++) {
		int row_max = maxloc(j, n, a);
		swap_row(n, a, b, j, row_max); 
		pthread_barrier_wait(&barrier_thread);
		block_size = (n - 1 - j) / p;
		start = (j + 1) + (block_size * (t_num - 1));
		for (i = start; i < min(start + block_size, n); i++) {
			double m = a[i][j] / a[j][j];
			for (k = j; k < n; k++) 
				a[i][k] -= m * a[j][k];
			b[i] -= m * b[j];
		}
		pthread_barrier_wait(&barrier_thread);
	}
	back_substitution(n, a, b, x);

	clock_gettime(CLOCK_REALTIME, &end);

	elapsed = SEC_TO_NANO(end.tv_sec - begin.tv_sec) + ((end.tv_nsec - begin.tv_nsec));
	printf("elapsed time: %lf (sec)\n", NANO_TO_SEC(elapsed));

	for (i = 0; i < n; i++) {
		double ax = 0;
		for (j = 0; j < n; j++)
			ax += a[i][j] * x[j];
		printf("%g\n", ax - b[i]);
	}

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

void *elimination(void *arg)
{
	struct args *aux = (struct args *) arg;
	int n = aux->n, p = aux->p;
	int t_num = aux->t_num;
	int step;
	int i, j;
	int block_size, start;
	double (*a)[n], (*b);

	a = (double (*)[n]) aux->a;
	b = aux->b;

	for (step = 0; step < n - 1; step++) {
		pthread_barrier_wait(&barrier_thread);
		block_size = (n - 1 - step) / p;
		start = (step + 1) + (block_size * (t_num - 1));
		/* Remained row is caculated by last thread */
		if (p == t_num) 
			block_size += (n - 1 - step) % p;
		for (i = start; i < min(start + block_size, n); i++) {
			double m = a[i][step] / a[step][step];
			for (j = step; j < n; j++) 
				a[i][j] -= m * a[step][j];
			b[i] -= m * b[step];
		}
		pthread_barrier_wait(&barrier_thread);	
	}
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
