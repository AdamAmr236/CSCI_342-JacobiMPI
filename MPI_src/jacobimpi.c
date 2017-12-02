/* Jacobi iteration using MPI

  usage ()

    mpicc -o jacobimpi jacobimpi.c 
    mpirun -np <numWorkers+1> jacobimpi  -- <gridSize> <iters>   */

#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#define MAXGRID 402     /* maximum grid size (real points plus edges) */
#define COORDINATOR 0   /* process number of the Coordinator */

static void Coordinator(int,int,int);
static void Worker(int,int,int,int,int);



/* main() -- initialize MPI, then become one of the worker or a coordinator  depending on your rank */

int main(int argc, char *argv[]) {
  int myid;
  int numWorkers, gridSize;  /* assume gridSize is multiple of numWorkers */
  int stripSize;             /* gridSize/numWorkers             */
  int numIters;              /* number of iterations to execute */
  double epsilon;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);  /* what is my id (rank)? */
  MPI_Comm_size(MPI_COMM_WORLD, &numWorkers);  /* how many processes? */
  numWorkers--;   /* one coordinator, the other processes are workers */

  /* get command-line arguments and do a simple error check */
  gridSize = atoi(argv[1]);
  //numIters = atoi(argv[2]);
  epsilon = atof(argv[2]);
  stripSize = gridSize/numWorkers;
  
  if (gridSize > (MAXGRID - 2) )
    {
      printf("Dumb sample program cannot handle a gridsize bigger than %d \n",MAXGRID);
      printf("Program does fixed size allocation to hold grid. \n");
      exit(1);
    }
  
  if (numWorkers < 2) {
    printf("Got to have at least two workers to have any fun.\n");
    printf("If there ended up beingonly one worker the dumb thing would still try to exchange \n");
    printf("messages with it's \"neighbors\" after each cycle\n");
    printf("Also serves as a crude check that we don't divide by zero computing the strip size. \n");
  }
  
  if (gridSize%numWorkers != 0) {
    printf("Dumb sample program cannot handle grid size that is not a multiple of the number of workers. \n");
    printf("Also note that the number of workers you give should actually be one more than a number that \n");
    printf("evenly divides the grid size as one of the workers is actually the coordinator. \n");
    printf("So if you want 4 actual workers to evnly divide a gridsize of 400, tell it 5 workers so the number\n");
    printf("of workers - 1 evenly divides the grid size. \n");
    printf("Also the gridsize you give does NOT include the borders. \n");
    exit(1);
  }
  /* Note it is an SPMD model so all the processes did the same thing up to this point now 
     become one of the workers or the coordinator, depending on my id */
  if (myid == 0) {
    printf("1 Coordinator and %d Workers\n", numWorkers);
    printf("  gridSize:  %d\n  stripSize:  %d\n  numIters:  %d\n",gridSize, stripSize, numIters);
    Coordinator(numWorkers, stripSize, gridSize);
  } else {
    Worker(myid, numWorkers, stripSize, gridSize, numIters);
  }

  MPI_Finalize();  /* clean up MPI */
}


/* gather and print results from Workers */
static void Coordinator(int numWorkers, int stripSize, int gridSize) {

  double grid[MAXGRID][MAXGRID];  // place to hold the results I get from the workers
  int i, j, startrow, endrow;
  int workerid;
  MPI_Status status;
  FILE *results;  // write the results to a file for viewing and debugging
  double mydiff = 0.0, maxdiff = 0.0;
  
  for (workerid = 1; workerid <= numWorkers; workerid++) {
    startrow = (workerid-1)*stripSize + 1;
    endrow = startrow + stripSize - 1;
    for (i = startrow; i <= endrow; i++) {
        MPI_Recv(&grid[i][1], gridSize, MPI_DOUBLE, workerid, 0,
            MPI_COMM_WORLD, &status);
    }
    printf("got results from worker %d\n", workerid);
  }
  // each work computes their maxdii and sends it to the coordinator to reduce
  MPI_Reduce(&mydiff, &maxdiff, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  printf("global maxdiff is %f\n", maxdiff);

  /* output the results to file "results" */
  results = fopen("results", "w");
  for (i = 1; i <= gridSize; i++) {
    for (j = 1; j <= gridSize; j++) {
      fprintf(results, "%f ", grid[i][j]);
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

static void Worker(int myid, int numWorkers, int stripSize,
                   int gridSize, int numIters) {
  /* the worker really only needs arrays of strip size + 2 not MAXGRID
     but since this dumb sample program does static allocation there is no way to know from 
     the command line parameter what strip size will end up being so it just wastes alot of space
     a "real" program (like the one you are writing) would dynamically allocate the memory 
  */
  double grid [2][MAXGRID][MAXGRID];
  int i, j, iters;
  int current = 0, next = 1;   /* current and next iteration indices */
  int left = 0, right = 0;    /* neighboring strips above and below */
  MPI_Status status;
  double mydiff = 0.0, maxdiff, temp;

  /* Your program should read the grid in from a file not just cram static values in like this 
     wimpy sample program 
  */
  
  /* set all points to 0s, then left and right edges of grid to 1s
     first worker sets top row to 1s; last worker sets bottom row to 1s */

  for (i = 0; i <= stripSize+1; i++)
    for (j = 0; j <= gridSize+1; j++) {
      grid[current][i][j] = 0.0;
      grid[next][i][j] = 0.0;
    }
  for (i = 0; i <= stripSize+1; i++) {
    grid[current][i][0] = 1.0;
    grid[current][i][gridSize+1] = 1.0;
    grid[next][i][0] = 1.0;
    grid[next][i][gridSize+1] = 1.0;
  }
  if (myid == 1)
    for (j = 0; j <= gridSize+1; j++) {
      grid[current][0][j] = 1.0;
      grid[next][0][j] = 1.0;
    }
  if (myid == numWorkers)
    for (j = 0; j <= gridSize+1; j++) {
      grid[current][stripSize+1][j] = 1.0;
      grid[next][stripSize+1][j] = 1.0;
    }

  /* determine neighbors */
  if (myid > 1)
    left = (myid-1);
  else
    left=numWorkers;
  
  if (myid < numWorkers)
    right = myid + 1;
  else
    right = 1;

  printf("Worker %d initialized; left is worker %d and right is worker %d\n", 
      myid, left, right);

  /* do the actual computation */
  for (iters = 1; iters <= numIters; iters++) {
    /* exchange my boundaries with my neighbors, in a ring */
    if (right != 0)
        MPI_Send(&grid[next][stripSize][1], gridSize, MPI_DOUBLE, right, 0,
                    MPI_COMM_WORLD);
    if (left != 0)
        MPI_Send(&grid[next][1][1], gridSize, MPI_DOUBLE, left, 0,
                    MPI_COMM_WORLD);
    if (left != 0)
        MPI_Recv(&grid[next][0][1], gridSize, MPI_DOUBLE, left, 0,
                    MPI_COMM_WORLD, &status);
    if (right != 0)
        MPI_Recv(&grid[next][stripSize+1][1], gridSize, MPI_DOUBLE, right, 0,
                    MPI_COMM_WORLD, &status);

    /* update my points */
    for (i = 1; i <= stripSize; i++) {
      for (j = 1; j <= gridSize; j++) {
        grid[next][i][j] = (grid[current][i-1][j] + grid[current][i+1][j] +
               grid[current][i][j-1] + grid[current][i][j+1]) / 4;
      }
    }

    /* swap roles of grids */
    current = next;  next = 1-next;
  }

  /* send results of my current strip to the coordinator */
  for (i = 1; i <= stripSize; i++) {
      MPI_Send(&grid[current][i][1], gridSize, MPI_DOUBLE,
            COORDINATOR, 0, MPI_COMM_WORLD);
  }

  /* compute maximum difference and reduce it with Coordinator */
  for (i = 1; i <= stripSize; i++)
    for (j = 1; j <= gridSize; j++) {
      temp = fabs(grid[next][i][j] - grid[current][i][j]);
      if (temp > mydiff)
        mydiff = temp;
    }
  MPI_Reduce(&mydiff, &maxdiff, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  printf("maxdiff of worker %d is %f\n", myid, mydiff);

}
