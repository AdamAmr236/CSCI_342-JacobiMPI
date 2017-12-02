#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#define  MASTER		0

int main (int argc, char *argv[])
{
unsigned int i,n;
 int rqd, pvd;
 rqd=MPI_THREAD_MULTIPLE;
 MPI_Init_thread(&argc, &argv,rqd,&pvd);
printf("Thread support level is %d ",pvd);

 switch (pvd)
   {
   case MPI_THREAD_SINGLE :
     printf("MPI_THREAD_SINGLE");
     break;
   case MPI_THREAD_FUNNELED :
     printf("MPI_THREAD_FUNNELED");
     break;
   case MPI_THREAD_SERIALIZED :
     printf("MPI_THREAD_SERIALIZED");
     break;
   case MPI_THREAD_MULTIPLE :
     printf("MPI_THREAD_MULTIPLE");
     break;
   }
 printf("\n");
MPI_Finalize();

}
