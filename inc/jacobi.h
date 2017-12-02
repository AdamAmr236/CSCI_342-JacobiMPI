#ifndef JACOBI
#define JACOBI

#define SHARED 1
#define BASE 10
#define MAXWORKERS 100

#define idx(x, i, j) ((x) * (i) + (j))

extern size_t x, y, xm1, ym1;
extern size_t firstY[MAXWORKERS];
extern size_t lastY[MAXWORKERS];
extern size_t numWorkers;
extern unsigned numIters;
extern double epsilon;
extern double maxDiff[MAXWORKERS];
extern double globalMaxDifference;
extern double * restrict grid1 __attribute__ ((aligned (16)));
extern double * restrict grid2 __attribute__ ((aligned (16)));

#endif
