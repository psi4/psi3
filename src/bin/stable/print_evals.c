/*! \file print_evals.c
    \ingroup (STABLE)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#define EXTERN
#include "globals.h"

void print_evals(double **evals, int *rank)
{
  int h, i;

  fprintf(outfile, "\t  #  ");
  for(h=0; h < moinfo.nirreps; h++)
      fprintf(outfile, "    %3s  ",moinfo.labels[h]);
  fprintf(outfile, "\n");
  fprintf(outfile, "\t---- ");
  for(h=0; h < moinfo.nirreps; h++)
      fprintf(outfile, "---------");
  fprintf(outfile, "\n");

  for(i=0; i < 5; i++) {
    fprintf(outfile, "\t %2d  ", i);
      for(h=0; h < moinfo.nirreps; h++) {
	  if(rank[h] <= i) fprintf(outfile, "         ");
	  else fprintf(outfile, " %7.4f ", evals[h][i]);
	}
      fprintf(outfile, "\n");
    }

  fprintf(outfile, "\n");
}
