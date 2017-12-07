/* Jacobi iteration using MPI

  usage ()

    mpicc -o jacobimpi jacobimpi.c
    mpirun -np <numWorkers+1> jacobimpi  -- <gridSize> <epsilon>   */

#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <assert.h>
#include <string.h>
#include <mcheck.h>
#define MAXGRID 402     /* maximum grid size (real points plus edges) */
#define COORDINATOR 0   /* process number of the Coordinator */
#define MAXTHREADS 8

#define idx(x, i, j) ((x) * (i) + (j))

static void Coordinator(int,int,int,double);
static void Worker(int,int,int,int,double);

double threadDiff[MAXTHREADS];
double globalMaxDiff;
double * restrict grid1 __attribute__ ((aligned (16)));
double * restrict grid2 __attribute__ ((aligned (16)));

int main(int argc, char *argv[]) {
  void *memptr __attribute__ ((aligned (16)));
  int mpiId, len;
  int numWorkers, gridSize;  /* assume gridSize is multiple of numWorkers */
  int stripSize;             /* gridSize/numWorkers             */
  int i, j;
  double epsilon;
  char hostname[MPI_MAX_PROCESSOR_NAME];

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiId);  /* what is my id (rank)? */
  MPI_Comm_size(MPI_COMM_WORLD, &numWorkers);  /* how many processes? */
  MPI_Get_processor_name(hostname, &len);
  numWorkers--;   /* one coordinator, the other processes are workers */
  // printf("Worker %d on %s\n", mpiId, hostname);

  /* get command-line arguments and do a simple error check */
  gridSize = atoi(argv[1]);
  size_t x = gridSize + 2;
  epsilon = atof(argv[2]);
  stripSize = gridSize/numWorkers;

  int n;
  if (mpiId != COORDINATOR){
    grid1 = malloc((x) * (stripSize+2) * sizeof(double));
    grid2 = malloc((x) * (stripSize+2) * sizeof(double));
  }else{
    grid1 = malloc((gridSize+2) * (gridSize+2) * sizeof(double));
  }

  if (numWorkers < 2) {
    printf("Must have at least two workers\n");
    exit(1);
  }

  if (gridSize%numWorkers != 0) {
    // printf("Dumb sample program cannot handle grid size that is not a multiple of the number of workers. \n");
    // printf("Also note that the number of workers you give should actually be one more than a number that \n");
    // printf("evenly divides the grid size as one of the workers is actually the coordinator. \n");
    // printf("So if you want 4 actual workers to evnly divide a gridsize of 400, tell it 5 workers so the number\n");
    // printf("of workers - 1 evenly divides the grid size. \n");
    // printf("Also the gridsize you give does NOT include the borders. \n");
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


  omp_set_num_threads(numWorkers);
  if (mpiId == COORDINATOR) {
    printf("1 Coordinator and %d Workers\n", numWorkers);
    printf("  gridSize:  %d\n  stripSize:  %d\n  epsilon:  %f\n",gridSize, stripSize, epsilon);
    Coordinator(numWorkers, stripSize, gridSize, epsilon);
  } else {
    #pragma omp parallel for
    for (int threadId = 0; threadId < 4; threadId++) { ///!!!!!!!!!!!!!!!!!!!!!!
      Worker(mpiId, numWorkers, stripSize, gridSize, epsilon);
    }
  }

  free(grid1);
  if (mpiId != COORDINATOR)
    free(grid2);
  printf("\nWorker %d Freed \n\n", mpiId);
  MPI_Finalize();  /* clean up MPI */
}


/* gather and print results from Workers */
static void Coordinator(int numWorkers, int stripSize, int gridSize, double epsilon) {

  int i, j, startrow, endrow;
  int workerid;
  MPI_Status status;
  FILE *results;  // write the results to a file for viewing and debugging
  double compDiff = 0.0;
  double maxdiff = epsilon + 1;
  int x = gridSize;

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
    printf("got results from worker %d\n", workerid);
  }
  // each work computes their maxdiff and sends it to the coordinator to reduce
  printf("global maxdiff is %f\n", maxdiff);

  /* output the results to file "results" */
  results = fopen("results", "w");
  for (i = 1; i <= gridSize; i++) {
    for (j = 1; j <= gridSize; j++) {
      fprintf(results, "%f ", grid1[(idx(x, i, j))]);
    }
    fprintf(results, "\n");
  }
}


/* Each Worker computes values in one strip, communicating with
  neighboring workers for border inforamation until the specified number of iterations.

  Your program, of course, runs till the overall compution converges which means that each
  worker can only know whether to stop or not based on the coordinator finding that everyones
  maxdiff is below the threshold.

  In the sample code we just get the maxdiff after each worker has done a fixed number of iterations
  Each then sends its results to the Coordinator for printing.
 */

static void Worker(int mpiId, int numWorkers, int stripSize,
                   int gridSize, double epsilon) {
  /* the worker really only needs arrays of strip size + 2 not MAXGRID
     but since this dumb sample program does static allocation there is no way to know from
     the command line parameter what strip size will end up being so it just wastes alot of space
     a "real" program (like the one you are writing) would dynamically allocate the memory
  */
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

  /* Your program should read the grid in from a file not just cram static values in like this
     wimpy sample program
  */

  /* set all points to 0s, then left and right edges of grid to 1s
     first worker sets top row to 1s; last worker sets bottom row to 1s */

  first = (threadId*(stripSize/numThreads)) + 1;
  last = (threadId+1) * (stripSize/numThreads);

  /* determine neighbors */
  if (mpiId > 1)
    left = (mpiId-1);
  else
    left=0;

  if (mpiId < numWorkers)
    right = mpiId + 1;
  else
    right = 0;

  #pragma omp single
  {
    printf("Worker %d initialized; left is worker %d and right is worker %d\n",
                mpiId, left, right);
  }

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
      if (right != 0)
          MPI_Send(&grid2[idx(x, stripSize, 1)], gridSize, MPI_DOUBLE, right, 0,
                      MPI_COMM_WORLD);
      if (left != 0)
          MPI_Send(&grid2[idx(x, 1, 1)], gridSize, MPI_DOUBLE, left, 0,
                      MPI_COMM_WORLD);
      if (left != 0)
          MPI_Recv(&grid2[idx(x, 0, 1)], gridSize, MPI_DOUBLE, left, 0,
                      MPI_COMM_WORLD, &status);
      if (right != 0)
          MPI_Recv(&grid2[idx(x, stripSize + 1, 1)], gridSize, MPI_DOUBLE, right, 0,
                      MPI_COMM_WORLD, &status);
    }

    #pragma omp barrier
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
  }

}
