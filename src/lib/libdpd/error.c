#include <stdio.h>

void dpd_error(char *caller, FILE *outfile)
{
  fprintf(outfile, "Error in: %s\n", caller);
  exit(2);
}
