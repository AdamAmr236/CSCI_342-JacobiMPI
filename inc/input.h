#ifndef JACOBI_INPUT
#define JACOBI_INPUT

void GetData(void);

unsigned long GetNumWorkers(const char *nptr);

unsigned long GetGridSize(const char *nptr);

double GetEpsilon(const char *nptr);

#endif
