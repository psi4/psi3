#include <stdio.h>

void idx_error(char *message, int p, int q, int r, int s, int pq, int rs,
	       int pq_sym, FILE *outfile)
{

  fprintf(outfile, "\n\tDPD Parameter Error in %s\n", message);
  fprintf(outfile,   "\t-------------------------------------------------\n");
  fprintf(outfile,   "\t    p      q      r      s  [   pq]  [   rs] symm\n");
  fprintf(outfile,   "\t%5d  %5d  %5d  %5d  [%5d]  [%5d]   %1d\n", p,q,r,s,
          pq,rs,pq_sym);
  exit(1);
}

