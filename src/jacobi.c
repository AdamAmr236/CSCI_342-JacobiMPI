/**
 * Jacobi iteration using pthreads
 */

#include <float.h>
#include <pthread.h>
#include <semaphore.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>

#include "barrier.h"
#include "input.h"
#include "jacobi.h"
#include "output.h"
#include "setup.h"
#include "worker.h"

size_t x;
size_t y;
size_t xm1;
size_t ym1;
size_t stripSize;
size_t numWorkers;
unsigned int numIters;
double epsilon;
double maxDiff[MAXWORKERS];
double *grid1;
double *grid2;

int main(int argc, char *argv[])
{
  pthread_t workerid[MAXWORKERS];
  pthread_attr_t attr;
  struct tms buffer;
  clock_t start, finish;
  size_t i;
  long threadId;
  double globalMaxDiff;

  if (argc != 5)
  {
    fprintf(stderr, "Wrong number of arguments.\n");
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "jacobi <x> <y> <numWorkers> <epsilon>\n");
    exit(EXIT_FAILURE);
  }

  x = GetGridSize(argv[1]);
  y = GetGridSize(argv[2]);
  numWorkers = GetNumWorkers(argv[3]);
  epsilon = GetEpsilon(argv[4]);

  xm1 = x - 1;
  ym1 = y - 1;
  stripSize = ym1 / numWorkers;
  numIters = 0;
  globalMaxDiff = DBL_MIN;
  grid1 = (double *) malloc(x * y * sizeof(double));
  grid2 = (double *) malloc(x * y * sizeof(double));

  InitializeGrids();

  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
  pthread_mutex_init(&barrier, NULL);
  pthread_cond_init(&go, NULL);

  start = times(&buffer);

#ifdef _DEBUG
  printf("Creating worker threads...\n");
#endif

  for (i = 0; i < numWorkers; i++)
  {
    threadId = i;
    pthread_create(&workerid[i], &attr, Worker, (void *) threadId);
  }

#ifdef _DEBUG
  printf("Waiting for all worker threads to finish...\n");
#endif

  for (i = 0; i < numWorkers; i++)
  {
    pthread_join(workerid[i], NULL);
  }

  finish = times(&buffer);

  for (i = 0; i < numWorkers; i++)
  {
    globalMaxDiff = (maxDiff[i] > globalMaxDiff) ? maxDiff[i] : globalMaxDiff;
  }

  printf("number of iterations:  %d\nmaximum difference:  %e\n", numIters, globalMaxDiff);
  printf("start:  %ld   finish:  %ld\n", start, finish);
  printf("elapsed time:  %ld\n", finish - start);

  WriteResults();

  exit(EXIT_SUCCESS);
}
