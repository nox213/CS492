#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>

#define NANO_TO_SEC(x) ((long) ((x) / 10e9))
#define SEC_TO_NANO(x) ((long) ((x) * 10e9))

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

int main(int argc, char *argv[])
{
	int n, p;
	int i, j;
	struct args aux;
	struct timespec begin, end;
	long elapsed;

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

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			a[i][j] = drand48();
			b[i][j] = drand48();
		}
	}

	init_task_queue(n);
	aux.a = a;
	aux.b = b;
	aux.c = c;
	aux.n = n;

	clock_gettime(CLOCK_REALTIME, &begin);
	for (i = 0; i < p; i++)
		pthread_create(&p_threads[i], NULL, mul_matrix, (void *) &aux);

	for (i = 0; i < p; i++)
		pthread_join(p_threads[i], NULL);
	clock_gettime(CLOCK_REALTIME, &end);
	
	elapsed = (end.tv_sec - begin.tv_sec) + (NANO_TO_SEC(end.tv_nsec - begin.tv_nsec));
	printf("elapsed time: %ld (sec)\n", (elapsed));

	return 0;
}

void print_array(int n, double arr[][n])
{
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) 
			printf("%lf ", arr[i][j]);
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
	
	a = aux->a;
	b = aux->b;
	c = aux->c;

	if ((local_sum = malloc(sizeof(double) * (n * n))) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return NULL;
	}
	memset(local_sum, 0, sizeof(double) * (n * n));

	while ((t = get_task()) != NULL) {
		int size = t->size;
		int sum;
		for (i = i_begin = t->i; i < i_begin + size; i++) {
			for (j = j_begin = t->j; j < j_begin + size; j++) {
				sum = 0;
				for (k = k_begin = t->j; k < k_begin + size; k++) {
					//printf("%d %d %d\n", i, j, k);
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
