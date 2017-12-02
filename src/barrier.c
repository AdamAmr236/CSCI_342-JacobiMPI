#include <pthread.h>
#include <stdio.h>
#include <float.h>

#include "barrier.h"
#include "jacobi.h"

pthread_mutex_t barrier;
pthread_cond_t go;
double globalMaxDifference;

int Barrier()
{
  static unsigned int numArrived = 0;

  pthread_mutex_lock(&barrier);

  numArrived++;

  if (numArrived == numWorkers)
  {
    numArrived = 0;

    /* increment iterations count */
    numIters++;
    
    /* swap grids */
    double *temp = grid1;
    grid1 = grid2;
    grid2 = temp; 

    /* find max of thread maxes */
    size_t i;
    globalMaxDifference = DBL_MIN;
    for (i = 0; i < numWorkers; i++) {
      globalMaxDifference = (maxDiff[i] > globalMaxDifference) ? maxDiff[i] : 
                                                        globalMaxDifference;
    }

    pthread_cond_broadcast(&go);
  }
  else
  {
    pthread_cond_wait(&go, &barrier);
  }

  pthread_mutex_unlock(&barrier);

  /* return 1 if convergence happened, 0 elsewise */
  if (globalMaxDifference <= epsilon) {
    return 1;
  } else {
    return 0;
  }
}
