#ifndef JACOBI
#define JACOBI

#define SHARED 1
#define BASE 10
#define MAXWORKERS 4

#define idx(x, i, j) ((x) * (i) + (j))

extern size_t x;
extern size_t y;
extern size_t xm1;
extern size_t ym1;
extern size_t stripSize;
extern size_t numWorkers;
extern unsigned int numIters;
extern double epsilon;
extern double maxDiff[MAXWORKERS];
extern double *grid1;
extern double *grid2;

#endif
