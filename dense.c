#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>

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


void print_array(int n, double arr[][n]);
void init_task_queue(int n);

int main(int argc, char *argv[])
{
	int n, p;
	int i, j, k;

	if (argc < 3) {
		fprintf(stderr, "dense n p\n");
		return 0;
	}

	n = strtol(argv[1], NULL, 10);
	p = strtol(argv[2], NULL, 10);

	double (*a)[n], (*b)[n], (*c)[n];
	pthread_t p_threads[n];

	if ((a = malloc(n * n)) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}
	if ((b = malloc(n * n)) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}
	if ((c = malloc(n * n)) == NULL) {
		fprintf(stderr, "malloc error: %d\n", __LINE__);
		return 0;
	}

	memset(p_threads, 0, sizeof(p_threads));
	memset(a, 0, n * n);
	memset(b, 0, n * n);
	memset(c, 0, n * n);

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			c[i][j] = drand48();
			b[i][j] = drand48();
		}
	}

	init_task_queue(n);

	/*
	print_array(n, b);
	print_array(n, c);
	*/

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

	pthread_mutex(&tq.queue_lock, NULL);
	tq.front = 0;
	if (n == 1000) {
		for (i = 0; i < 10; i++)  {
			for (j  = 0; j < 10; j++) {
				tq.t[i].size = 100;
				tq.t[i].i = i * 100;
				tq.t[i].j = j * 100;
			}
		}
	}
	else if (n == 3000) {
		for (i = 0; i < 10; i++)  {
			for (j  = 0; j < 10; j++) {
				tq.t[i].size = 300;
				tq.t[i].i = i * 300;
				tq.t[i].j = j * 300;
			}
		}
	}
	else {
		for (i = 0; i < 10; i++)  {
			for (j  = 0; j < 10; j++) {
				tq.t[i].size = 500;
				tq.t[i].i = i * 500;
				tq.t[i].j = j * 500;
			}
		}
	}
}

struct task *get_task(void)
{
	struct task *tp = NULL;

	pthread_mutex_lock(&tq.queue_lock);
	if (front < 100)
		tp = &tq.t[front++];
	pthread_mutex_unlock(&tq.queue_lock);

	return tp;
}

void *mul_matrix(void *arg)
{

}



