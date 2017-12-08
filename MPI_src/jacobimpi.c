/* Jacobi iteration using MPI

  usage ()

    mpicc -fopenmp -o jacobimpi jacobimpi.c
    mpirun -bynode -np <numWorkers+1> --hostfile <hostfile> jacobimpi <gridSize> <numThreads> <epsilon>   */

#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <time.h>
#define COORDINATOR 0   /* process number of the Coordinator */
#define MAXTHREADS 100

#define idx(x, i, j) ((x) * (i) + (j))

static void Coordinator(int,int,int,double);
static void Worker(int,int,int,int,double);

double threadDiff[MAXTHREADS];
double globalMaxDiff;
double * grid1;
double * grid2;

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

int main(int argc, char *argv[]) {
  int mpiId, numThreads, len;
  int numWorkers, gridSize;  /* assume gridSize is multiple of numWorkers */
  int stripSize;             /* gridSize/numWorkers             */
  int i, j;
  double epsilon;
  char hostname[MPI_MAX_PROCESSOR_NAME];
  FILE *results;  // write the results to a file for viewing and debugging
  clock_t start, finish;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiId);  /* what is my id (rank)? */
  MPI_Comm_size(MPI_COMM_WORLD, &numWorkers);  /* how many processes? */
  MPI_Get_processor_name(hostname, &len);
  numWorkers--;   /* one coordinator, the other processes are workers */
  printf("Worker %d on %s\n", mpiId, hostname);
  MPI_Barrier(MPI_COMM_WORLD);



  /* get command-line arguments and do a simple error check */
  gridSize = atoi(argv[1]);
  size_t x = gridSize + 2;
  numThreads = atoi(argv[2]);
  epsilon = atof(argv[3]);
  stripSize = gridSize/numWorkers;

  if (mpiId == COORDINATOR)
    clock_gettime(CLOCK_REALTIME, &clockStart);

  int n;
  if (mpiId != COORDINATOR){
    if (mpiId == numWorkers) {
      /* edge case for when gridSize isn't evenly divisible by numWorkers*/
      grid1 = malloc((x) * (stripSize + 2 + (gridSize % numWorkers)) * sizeof(double));
      grid2 = malloc((x) * (stripSize + 2 + (gridSize % numWorkers)) * sizeof(double));
    } else {
      grid1 = malloc((x) * (stripSize+2) * sizeof(double));
      grid2 = malloc((x) * (stripSize+2) * sizeof(double));
    }
  } else {
    grid1 = malloc((gridSize+2) * (gridSize+2) * sizeof(double));
  }

  if (numWorkers < 2) {
    printf("Must have at least two workers\n");
    exit(1);
  }


  /* Note it is an SPMD model so all the processes did the same thing up to this point now
     become one of the workers or the coordinator, depending on my id */

  if (mpiId != COORDINATOR){
    int current = 0, next = 1;
    for (i = 0; i <= stripSize+1; i++)
     for (j = 0; j <= gridSize+1; j++) {
       grid1[idx(x, i, j)] = 0.0;
       grid2[idx(x, i, j)] = 0.0;
     }
    for (i = 0; i <= stripSize+1; i++) {
     grid1[idx(x, i, 0)] = 1.0;
     grid1[idx(x, i, gridSize+1)] = 1.0;
     grid2[idx(x, i, 0)] = 1.0;
     grid2[idx(x, i, gridSize+1)] = 1.0;
    }
    if (mpiId == 1)
     for (j = 0; j <= gridSize+1; j++) {
       grid1[idx(x, 0, j)] = 1.0;
       grid2[idx(x, 0, j)] = 1.0;
     }
    if (mpiId == numWorkers)
     for (j = 0; j <= gridSize+1; j++) {
       grid1[idx(x, stripSize + 1, j)] = 1.0;
       grid2[idx(x, stripSize + 1, j)] = 1.0;
     }
   }

   //printf("omp_get_max_num_threads: %d\n", omp_get_max_threads());
  omp_set_num_threads(numThreads);
  if (mpiId == COORDINATOR) {
    printf("1 Coordinator and %d Workers\n", numWorkers);
    printf("  gridSize:  %d\n  stripSize:  %d\n  epsilon:  %f\n",gridSize, stripSize, epsilon);
    Coordinator(numWorkers, stripSize, gridSize, epsilon);
  } else {
    #pragma omp parallel
      { Worker(mpiId, numWorkers, stripSize, gridSize, epsilon); }
  }

  if (mpiId == COORDINATOR)
    clock_gettime(CLOCK_REALTIME, &clockFinish);

  struct timespec clockTime = diff(&clockStart, &clockFinish);
  if (mpiId == COORDINATOR) {
    printf("Wall Clock Time: %ld.%ld\n\n",
              clockTime.tv_sec, clockTime.tv_nsec);
  }

  /* output the results to file "results" */
  if (mpiId == COORDINATOR) {
    results = fopen("results", "w");
    for (i = 1; i <= gridSize; i++) {
      for (j = 1; j <= gridSize; j++) {
        fprintf(results, "%f ", grid1[(idx(x, i, j))]);
      }
      fprintf(results, "\n");
    }
  }
  free(grid1);
  if (mpiId != COORDINATOR)
  free(grid2);


  MPI_Finalize();  /* clean up MPI */
}


/* gather and print results from Workers */
static void Coordinator(int numWorkers, int stripSize, int gridSize, double epsilon) {

  int i, j, startrow, endrow;
  int workerid;
  MPI_Status status;
  double compDiff = 0.0;
  double maxdiff = epsilon + 1;
  size_t x = gridSize + 2;

  while (maxdiff > epsilon) {
    MPI_Allreduce(&compDiff, &maxdiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }

  for (workerid = 1; workerid <= numWorkers; workerid++) {
    startrow = (workerid-1)*stripSize + 1;
    endrow = startrow + stripSize - 1;
    for (i = startrow; i <= endrow; i++) {
        MPI_Recv(&grid1[idx(x, i, 1)], gridSize, MPI_DOUBLE, workerid, 0,
            MPI_COMM_WORLD, &status);
    }
  }
  // each work computes their maxdiff and sends it to the coordinator to reduce
  printf("global maxdiff is %f\n", maxdiff);

}

static void Worker(int mpiId, int numWorkers, int stripSize,
                   int gridSize, double epsilon) {
  int i, j;
  int threadId = omp_get_thread_num();
  int numThreads = omp_get_num_threads();
  int current = 0, next = 1;   /* current and next iteration indices */
  int left = 0, right = 0;     /* neighboring strips above and below */
  int first, last;             /* first and last strip for omp */
  MPI_Status status;
  double maxdiff, temp;
  double diff;
  int x = gridSize + 2;


  first = (threadId*(stripSize/numThreads)) + 1;
  if (threadId == numThreads - 1) {
    last = stripSize;
  } else {
    last = (threadId+1) * (stripSize/numThreads);
  }

  /* determine neighbors */
  if (mpiId > 1)
    left = (mpiId-1);
  else
    left=0;

  if (mpiId < numWorkers)
    right = mpiId + 1;
  else
    right = 0;


  // Each thread has their own maxdiff
  // Each computer has one globalMaxDiff
  maxdiff = epsilon + 1;
  #pragma omp single
  {
    globalMaxDiff = epsilon + 1;
  }

  while (globalMaxDiff > epsilon) {

    // Calculate Jacobi iteration
    threadDiff[threadId] = 0;
    for (i = first; i <= last; i++) {
      for (j = 1; j <= gridSize; j++) {
        grid2[idx(x, i, j)] = (grid1[idx(x, i-1, j)] + grid1[idx(x, i + 1, j)] +
               grid1[idx(x, i, j-1)] + grid1[idx(x, i, j+1)]) * 0.25;
        diff = grid2[idx(x, i, j)] - grid1[idx(x, i, j)];
        if (diff < 0) {
          diff = - diff;
        }
        if (diff > threadDiff[threadId]) {
          threadDiff[threadId] = diff;
        }
      }
    }

    // Send borders to each other
    #pragma omp barrier
    #pragma omp single
    {
      if (right != 0) {
        MPI_Send(&grid2[idx(x, stripSize, 1)], gridSize, MPI_DOUBLE, right, 0,
                  MPI_COMM_WORLD);
      }

      if (left != 0) {
        MPI_Recv(&grid2[idx(x, 0, 1)], gridSize, MPI_DOUBLE, left, 0,
                  MPI_COMM_WORLD, &status);
      }

      if (left != 0) {
        MPI_Send(&grid2[idx(x, 1, 1)], gridSize, MPI_DOUBLE, left, 0,
                  MPI_COMM_WORLD);
      }

      if (right != 0) {
        MPI_Recv(&grid2[idx(x, stripSize + 1, 1)], gridSize, MPI_DOUBLE, right, 0,
                  MPI_COMM_WORLD, &status);
      }
    }

    //printf("Worker %d - Thread %d at barrier\n", mpiId, threadId);
    #pragma omp single
    {
      // Swap Grids
      double* temp = grid1;
      grid1 = grid2;
      grid2 = temp;

      // Calculate maxdiff for each computer
      globalMaxDiff = -1;
      size_t k;
      for (k = 0; k < numWorkers; k++){
        globalMaxDiff = (threadDiff[k] > globalMaxDiff) ?
          threadDiff[k] : globalMaxDiff;
      }
      // Calculate maxdiff between each computer, set to
      // each computers globalMaxDiff
      MPI_Allreduce(&globalMaxDiff, &maxdiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      globalMaxDiff = maxdiff;
    }
  }

  #pragma omp barrier
  #pragma omp single
  {
    for (i = 1; i <= stripSize; i++) {
        MPI_Send(&grid1[idx(x, i, 1)], gridSize, MPI_DOUBLE,
              COORDINATOR, 0, MPI_COMM_WORLD);
    }
    //printf("finished worker thread\n");
  }

}
