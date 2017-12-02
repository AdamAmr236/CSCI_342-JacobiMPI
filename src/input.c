#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "input.h"
#include "jacobi.h"

static unsigned long GetUnsignedLong(const char *nptr)
{
  unsigned long result;

  errno = 0;

  result = strtoul(nptr, NULL, BASE);

  if (errno != 0)
  {
    fprintf(stderr, "%s\n", strerror(errno));
    exit(EXIT_FAILURE);
  }

  return result;
}

static double GetDouble(const char *nptr)
{
  double result;

  errno = 0;

  result = strtod(nptr, NULL);

  if (errno != 0)
  {
    fprintf(stderr, "%s\n", strerror(errno));
    exit(EXIT_FAILURE);
  }

  return result;
}

void GetData()
{
  char *line = NULL;
  size_t length = 0;
  ssize_t count;
  const char *delim = " ";
  char *nptr, *endptr;
  size_t i, j;
  double current;

  for (i = 0; i < y; i++)
  {
    count = getline(&line, &length, stdin);

    if (count == EOF && ferror(stdin) != 0)
    {
      fprintf(stderr, "Failed to read input.\n");
      exit(EXIT_FAILURE);
    }

    nptr = strtok(line, delim);

    for (j = 0; nptr != NULL; j++)
    {
      current = strtod(nptr, &endptr);

      if (nptr != endptr)
      {
        grid1[idx(x, i, j)] = grid2[idx(x, i, j)] = current;
      }

      nptr = strtok(NULL, delim);
    }
  }

  free(line);
}

unsigned long GetNumWorkers(const char *nptr)
{
  unsigned long result = GetUnsignedLong(nptr);

  if (result < 1)
  {
    fprintf(stderr, "Thread count cannot be less than 1.\n");
    exit(EXIT_FAILURE);
  }

  return result;
}

unsigned long GetGridSize(const char *nptr)
{
  unsigned long result = GetUnsignedLong(nptr);

  if (result < 2)
  {
    fprintf(stderr, "Grid size cannot be less than 1.\n");
    exit(EXIT_FAILURE);
  }

  return result;
}

double GetEpsilon(const char *nptr)
{
  double result = GetDouble(nptr);

  if (result <= 0)
  {
    fprintf(stderr, "Epsilon must be greater than 0.\n");
    exit(EXIT_FAILURE);
  }

  return result;
}
