#include <stdlib.h>

#include "jacobi.h"
#include "setup.h"

void InitializeGrids()
{
  size_t i, j;

  for (i = 0; i < y; i++)
  {
    for (j = 0; j < x; j++)
    {
      grid1[idx(x, i ,j)] = 0.0;
      grid2[idx(x, i, j)] = 0.0;
    }
  }

  for (i = 0; i < y; i++)
  {
    grid1[idx(x, i, 0)] = 1.0;
    grid1[idx(x, i, xm1)] = 1.0;
    grid2[idx(x, i, 0)] = 1.0;
    grid2[idx(x, i, xm1)] = 1.0;
  }

  for (j = 0; j < x; j++)
  {
    grid1[idx(x, 0, j)] = 1.0;
    grid2[idx(x, 0, j)] = 1.0;
    grid1[idx(x, ym1, j)] = 1.0;
    grid2[idx(x, ym1, j)] = 1.0;
  }
}
