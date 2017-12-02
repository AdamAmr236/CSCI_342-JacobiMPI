#include <stdio.h>

#include "jacobi.h"
#include "output.h"

void WriteResults(FILE *results)
{
  size_t i, j;

  for (i = 1; i < ym1; i++)
  {
    for (j = 1; j < xm1; j++)
    {
      fprintf(results, "%f ", grid2[idx(x, i, j)]);
    }

    fprintf(results, "\n");
  }
}
