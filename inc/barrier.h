#ifndef JACOBI_BARRIER
#define JACOBI_BARRIER

extern pthread_mutex_t barrier;
extern pthread_cond_t go;

int Barrier(void);

#endif
