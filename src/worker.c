#include <float.h>
#include <pthread.h>
#include <stdio.h>

#include "barrier.h"
#include "jacobi.h"
#include "worker.h"

void *Worker(void *arg)
{
  long id = (long) arg;
  size_t first, last, i, j;
  int iters;
  double diff, localMaxDiff, globalMaxDiff;

#ifdef _DEBUG
  printf("worker %ld pthread id %lu has started\n", id, pthread_self());
#endif

  first = id * stripSize + 1;
  last = first + stripSize - 1;
  globalMaxDiff = localMaxDiff = DBL_MAX;

  while (globalMaxDiff > epsilon)
  {
    for (i = first; i < last; i++)
    {
      for (j = 1; j < xm1; j++)
      {
        grid2[idx(x, i, j)] = (grid1[idx(x, i - 1, j)] + grid1[idx(x, i + 1, j)] + grid1[idx(x, i, j - 1)] + grid1[idx(x, i, j + 1)]) * 0.25;
      }
    }

    iters++;

    Barrier();

    globalMaxDiff = localMaxDiff = DBL_MIN;

    for (i = first; i < last; i++)
    {
      for (j = 1; j < xm1; j++)
      {
        grid1[idx(x, i, j)] = (grid2[idx(x, i - 1, j)] + grid2[idx(x, i + 1, j)] + grid2[idx(x, i, j - 1)] + grid2[idx(x, i, j + 1)]) * 0.25;
        diff = grid1[idx(x, i, j)] - grid2[idx(x, i, j)];
        diff = (diff < 0) ? -diff : diff;
        localMaxDiff = (diff > localMaxDiff) ? diff : localMaxDiff;
      }
    }

    maxDiff[id] = localMaxDiff;

    iters++;

    Barrier();

    for (i = 0; i < numWorkers; i++)
    {
      globalMaxDiff = (maxDiff[i] > globalMaxDiff) ? maxDiff[i] : globalMaxDiff;
    }
  }

#ifdef _DEBUG
  printf("worker %ld pthread id %lu has finished\n", id, pthread_self());
#endif

  return NULL;
}
