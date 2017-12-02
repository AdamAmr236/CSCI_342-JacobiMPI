#include <stdio.h>

#include "jacobi.h"
#include "output.h"

void WriteResults()
{
  size_t i, j;
  FILE *results = fopen("results", "w");

  for (i = 1; i < y; i++)
  {
    for (j = 1; j < x; j++)
    {
      fprintf(results, "%f ", grid2[idx(x, i, j)]);
    }

    fprintf(results, "\n");
  }
}
