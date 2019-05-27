#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <stdint.h>

#define SEC_TO_NANO(x) (((x) * 1000000000L))

struct task {
	int i;
	int j;
	int size;
};

struct task_queue {
	struct task t[100];
	int front;
	pthread_mutex_t queue_lock;
} tq;

struct args {
	double **a, **b, **c;
	int n;
};

pthread_mutex_t c_lock = PTHREAD_MUTEX_INITIALIZER;

void print_array(int n, double arr[][n]);
void init_task_queue(int n);
void *mul_matrix(void *arg);
void mul_matrix_single(int n, double (*a)[n], double (*b)[n], double (*answer)[n]);
bool nearly_equal(double a, double b, double epsilon); 

static inline double min(const double a, const double b)
{
	return a < b ? a : b;
}

int main(int argc, char *argv[])
{
	int n, p;
	int i, j;
	struct args aux;
	struct timespec begin, end;
	uint64_t elapsed_s, elapsed_p;

	if (argc < 3) {
		fprintf(stderr, "dense n p\n");
		return 0;
	}

	n = strtol(argv[1], NULL, 10);
	p = strtol(argv[2], NULL, 10);

	double (*a)[n], (*b)[n], (*c)[n], (*answer)[n];
	pthread_t p_threads[p];

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

	memset(p_threads, 0, sizeof(p_threads));
	memset(a, 0, sizeof(double) * (n * n));
	memset(b, 0, sizeof(double) * (n * n));
	memset(c, 0, sizeof(double) * (n * n));
	memset(answer, 0, sizeof(double) * (n * n));

	srand48(time(NULL));
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			a[i][j] = drand48();
			b[i][j] = drand48();
		}
	}


	printf("Multi thread computaion start\n");
	clock_gettime(CLOCK_REALTIME, &begin);

	init_task_queue(n);
	aux.a = (double **) a;
	aux.b = (double **) b;
	aux.c = (double **) c;
	aux.n = n;

	if (p == 1)
		mul_matrix_single(n, a, b, c);
	else  {

		for (i = 0; i < p; i++)
			pthread_create(&p_threads[i], NULL, mul_matrix, (void *) &aux);

		for (i = 0; i < p; i++)
			pthread_join(p_threads[i], NULL);
	}
	clock_gettime(CLOCK_REALTIME, &end);

	elapsed_p = SEC_TO_NANO(end.tv_sec - begin.tv_sec) + (end.tv_nsec - begin.tv_nsec);
	printf("Multi thread computaion end\n");
	printf("elapsed time: %ld (nsec)\n", elapsed_p);

	printf("Single thread computaion start\n");
	clock_gettime(CLOCK_REALTIME, &begin);
	mul_matrix_single(n, a, b, answer);
	clock_gettime(CLOCK_REALTIME, &end);
	elapsed_s = SEC_TO_NANO(end.tv_sec - begin.tv_sec) + (end.tv_nsec - begin.tv_nsec);
	printf("Single thread computaion end\n");
	printf("elapsed time: %ld (nsec)\n", elapsed_s);

	printf("speed up is %g\n", (double) elapsed_s / elapsed_p);

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

void print_array(int n, double arr[][n])
{
	int i, j;

	for (i = 0; i < 1; i++) {
		for (j = 0; j < n; j++) 
			printf("%g ", arr[i][j]);
		printf("\n");
	}
}

void init_task_queue(int n)
{
	int i, j;

	pthread_mutex_init(&tq.queue_lock, NULL);
	tq.front = 0;
	if (n == 1000) {
		for (i = 0; i < 10; i++)  {
			for (j  = 0; j < 10; j++) {
				tq.t[i*10+j].size = 100;
				tq.t[i*10+j].i = i * 100;
				tq.t[i*10+j].j = j * 100;
			}
		}
	}
	else if (n == 3000) {
		for (i = 0; i < 10; i++)  {
			for (j  = 0; j < 10; j++) {
				tq.t[i*10+j].size = 300;
				tq.t[i*10+j].i = i * 300;
				tq.t[i*10+j].j = j * 300;
			}
		}
	}
	else {
		for (i = 0; i < 10; i++)  {
			for (j  = 0; j < 10; j++) {
				tq.t[i*10+j].size = 500;
				tq.t[i*10+j].i = i * 500;
				tq.t[i*10+j].j = j * 500;
			}
		}
	}
}

struct task *get_task(void)
{
	struct task *tp = NULL;

	pthread_mutex_lock(&tq.queue_lock);
	if (tq.front < 100)
		tp = &tq.t[tq.front++];
	pthread_mutex_unlock(&tq.queue_lock);

	return tp;
}

void *mul_matrix(void *arg)
{
	int i, j, k;
	int i_begin, j_begin, k_begin;
	struct args *aux = (struct args *) arg;
	int n = aux->n;
	struct task *t;
	double (*local_sum)[n], (*a)[n], (*b)[n], (*c)[n];

	a = (double (*)[n]) aux->a;
	b = (double (*)[n]) aux->b;
	c = (double (*)[n]) aux->c;

	if ((local_sum = malloc(sizeof(double) * (n * n))) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return NULL;
	}
	memset(local_sum, 0, sizeof(double) * (n * n));

	while ((t = get_task()) != NULL) {
		int size = t->size;
		for (i = 0; i < n; i++) {
			for (j = j_begin = t->i; j < j_begin + size; j++) {
				double sum = 0;
				for (k = k_begin = t->j; k < k_begin + size; k++) {
					sum += a[i][k] * b[k][j];
				}
				local_sum[i][j] += sum;
			}
		}
	}


	pthread_mutex_lock(&c_lock);
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			c[i][j] += local_sum[i][j];
	pthread_mutex_unlock(&c_lock);

	return NULL;
}

void mul_matrix_single(int n, double (*a)[n], double (*b)[n], double (*answer)[n])
{
	int i, j, k;

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			for (k = 0; k < n; k++)
				answer[i][j] += a[i][k] * b[k][j];
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
