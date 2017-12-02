#include <pthread.h>
#include <stdio.h>

#include "barrier.h"
#include "jacobi.h"

pthread_mutex_t barrier;
pthread_cond_t go;

void Barrier()
{
  static unsigned int numArrived = 0;

  pthread_mutex_lock(&barrier);

  numArrived++;

  if (numArrived == numWorkers)
  {
    numArrived = 0;
    pthread_cond_broadcast(&go);
  }
  else
  {
    pthread_cond_wait(&go, &barrier);
  }

  pthread_mutex_unlock(&barrier);
}
