/**
 * Jacobi iteration using openMP
 */

#include <float.h>
#include <math.h>
#include <semaphore.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/times.h>
#include <time.h>
#include <omp.h>

#include "barrier.h"
#include "input.h"
#include "jacobi.h"
#include "output.h"
#include "setup.h"
#include "worker.h"
size_t x, y, xm1, ym1;
size_t firstY[MAXWORKERS];
size_t lastY[MAXWORKERS];
size_t numWorkers;
unsigned numIters;
double epsilon;
double maxDiff[MAXWORKERS];
double * restrict grid1 __attribute__ ((aligned (16)));
double * restrict grid2 __attribute__ ((aligned (16)));
int squaresize;
struct List *squareLists;


void setupSquareSlices(void);

/* used for timing */
struct timespec clockStart, clockFinish;

struct timespec diff(struct timespec *startTime, struct timespec *endTime)
{
  struct timespec temp;
  if ((endTime->tv_nsec-startTime->tv_nsec)<0) {
    temp.tv_sec = endTime->tv_sec - startTime->tv_sec-1;
    temp.tv_nsec = 1000000000 + endTime->tv_nsec - startTime->tv_nsec;
  } else {
    temp.tv_sec = endTime->tv_sec - startTime->tv_sec;
    temp.tv_nsec = endTime->tv_nsec - startTime->tv_nsec;
  }
  return temp;
}

int main(int argc, char *argv[])
{
  void *memptr __attribute__ ((aligned (16)));
  int status;
  struct tms buffer;
  clock_t start, finish;
  size_t i;
  uint64_t cyclesHighStart, cyclesLowStart, cyclesHighEnd, cyclesLowEnd, nstart, nfinish;
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
  numIters = 0;
  globalMaxDiff = DBL_MIN;

  int alignment = 16;

  status = posix_memalign(&memptr, alignment, x * y * sizeof(double));

  if (status != 0)
  {
    fprintf(stderr, "%s\n", strerror(status));
    exit(EXIT_FAILURE);
  }
  else
  {
    grid1 = (double *) __builtin_assume_aligned(memptr, alignment);
    //grid1 = calloc(x * y, sizeof(double));
  }

  status = posix_memalign(&memptr, alignment, x * y * sizeof(double));

  if (status != 0)
  {
    fprintf(stderr, "%s\n", strerror(status));
    exit(EXIT_FAILURE);
  }
  else
  {
    grid2 = (double *) __builtin_assume_aligned(memptr, alignment);
    //grid2 = calloc(x * y, sizeof(double));
  }

  InitializeGrids();
  CalculateSlices();

  clock_gettime(CLOCK_REALTIME, &clockStart);

  asm volatile (
      "RDTSC\n\t"
      "movq %%rdx, %0\n\t"
      "movq %%rax, %1\n\t": "=r" (cyclesHighStart), "=r" (cyclesLowStart)::
    "%rax", "%rbx", "%rcx", "%rdx");

  start = times(&buffer);

#ifdef _DEBUG
  printf("Creating worker threads...\n");
#endif

  omp_set_num_threads(numWorkers);
  #pragma omp parallel for
  for (i = 0; i < numWorkers; i++)
  {
    Worker(i);
  }

#ifdef _DEBUG
  printf("Waiting for all worker threads to finish...\n");
#endif

clock_gettime(CLOCK_REALTIME, &clockFinish);

asm volatile (
    "RDTSC\n\t"
    "movq %%rdx, %0\n\t"
    "movq %%rax, %1\n\t": "=r" (cyclesHighEnd), "=r" (cyclesLowEnd)::
  "%rax", "%rbx", "%rcx", "%rdx");

  nstart = (((uint64_t) cyclesHighStart << 32) | cyclesLowStart);
  nfinish = (((uint64_t) cyclesHighEnd << 32) | cyclesLowEnd);

  finish = times(&buffer);

  for (i = 0; i < numWorkers; i++)
  {
    globalMaxDiff = fmax(globalMaxDiff, maxDiff[i]);
  }

  // printf("number of iterations:  %d\nmaximum difference:  %e\n", numIters, globalMaxDiff);
  // printf("start:  %ld   finish:  %ld\n", start, finish);
  // printf("elapsed time:  %ld\n", finish - start);
  printf("start:  %llu   end:  %llu\n", (long long unsigned) nstart, (long long unsigned) nfinish);
  printf("number of cycles:  %llu\n", (long long unsigned) (nfinish - nstart));

  /* get difference */
  struct timespec clockTime = diff(&clockStart, &clockFinish);
  printf("wall clock time: %ld.%ld\n", clockTime.tv_sec, clockTime.tv_nsec);
  printf("\n");

  //printf("Middle of array:  %f\n", grid1[idx(x, 500, 500)]);

  //results = fopen("results", "w");

  //WriteResults(stdout);

  exit(EXIT_SUCCESS);
}
