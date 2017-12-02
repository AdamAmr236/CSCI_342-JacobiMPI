#ifndef JACOBI_WORKER
#define JACOBI_WORKER

void Worker(size_t i);

void *WorkerSquares(void *);


struct Coords {
    size_t x;
    size_t y;
};

struct List {
    struct Coords *array;
    size_t size;
};


extern int squaresize;
extern struct List *squareLists;
#endif
