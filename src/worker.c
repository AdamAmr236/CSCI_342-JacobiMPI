#include <float.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <omp.h>

#include "barrier.h"
#include "jacobi.h"
#include "worker.h"

void Worker(size_t id)
{
  size_t first, last, i, j;
  double localMaxDiff;

#ifdef _DEBUG
  printf("worker %lu pthread self_id %lu has started\n", id, pthread_self());
#endif

  first = firstY[id];
  last = lastY[id];
  #pragma omp single
  {
    globalMaxDifference = localMaxDiff = DBL_MAX; 
  }

  while (globalMaxDifference > epsilon)
  {
    localMaxDiff = DBL_MIN;

    for (i = first; i < last; i++)
    {
      for (j = 1; j < xm1; j++)
      {
        grid1[idx(x, i, j)] = (grid2[idx(x, i - 1, j)] + grid2[idx(x, i + 1, j)] + grid2[idx(x, i, j - 1)] + grid2[idx(x, i, j + 1)]) * 0.25;
        localMaxDiff = fmax(localMaxDiff, fabs(grid1[idx(x, i, j)] - grid2[idx(x, i, j)]));
      }
    }

    maxDiff[id] = localMaxDiff;

    #pragma omp barrier
    #pragma omp single
    {
      /* increment iterations count */
      numIters++;
      //printf("Iter: %u", numIters);
      
      /* swap grids */
      double *temp = grid1;
      grid1 = grid2;
      grid2 = temp; 

      /* find max of thread maxes */
      size_t k;
      globalMaxDifference = DBL_MIN;
      for (k = 0; k < numWorkers; k++) {
        //printf("%f\n", maxDiff[k]);
        globalMaxDifference = (maxDiff[k] > globalMaxDifference) ? maxDiff[k] : 
                                                          globalMaxDifference;
      }
    }
  }

#ifdef _DEBUG
  printf("worker %ld pthread id %lu has finished\n", id, pthread_self());
#endif
}

void *WorkerSquares(void *arg)
{
  long id = (long) arg;
  int converged = 0;
  double diff;


  #ifdef _DEBUG
    printf("worker %ld pthread id %lu has started\n", id, pthread_self());
  #endif

  struct List list = squareLists[id];

  //printf("Coords: %lu, %lu\n", list.array[2].x, list.array[2].y);

  while (!converged) {
    maxDiff[id] = DBL_MIN;
    /* for every (x,y) coordinate pair in the list */
    for (size_t k=0; k < list.size; k++) {
      //printf("K=%lu\n", k);
      struct Coords coords = list.array[k];

      size_t rowCutoff = (coords.y + squaresize) < ym1 ? (coords.y + squaresize) : ym1; 
      size_t colCutoff = (coords.x + squaresize) < xm1 ? (coords.x + squaresize) : xm1; 

      for (size_t i = coords.y; i < rowCutoff; i++) {
        for (size_t j = coords.x; j < colCutoff; j++) {
          grid2[idx(x, i, j)] = (grid1[idx(x, i - 1, j)] + grid1[idx(x, i + 1, j)] + grid1[idx(x, i, j - 1)] + grid1[idx(x, i, j + 1)]) * 0.25;
          /* check difference and save max in global array*/
          diff = grid2[idx(x, i, j)] - grid1[idx(x, i, j)];
          diff = (diff < 0) ? -diff : diff;
          maxDiff[id] = (diff > maxDiff[id]) ? diff : maxDiff[id];
        }
      }
    }
  
    /* the barrier will check for convergance after all threads check in */
    converged = Barrier();
  
  }


  #ifdef _DEBUG
    printf("worker %ld pthread id %lu has finished\n", id, pthread_self());
  #endif

  return NULL;
}