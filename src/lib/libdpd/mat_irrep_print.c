#include <stdio.h>
#include <stdlib.h>
#include "dpd.h"

int dpd_mat_irrep_print(double **matrix, struct dpdparams *Params,
			int irrep, FILE *outfile)
{
  div_t fraction;
  int i,j;
  int rows, cols, cols_per_page, num_pages, last_page, page, first_col;

  cols_per_page = 5;

  rows = Params->rowtot[irrep];
  cols = Params->coltot[irrep];

  /* Determine the number of cols_per_page groups */
  fraction = div(cols,cols_per_page);
  num_pages = fraction.quot;  /* Number of complete column groups */
  last_page = fraction.rem;  /* Number of columns in last group */

  /* Loop over the complete column groups */
  for(page=0; page < num_pages; page++) {
      first_col = page*cols_per_page;

      fprintf(outfile,"\n           ");
      for(i=first_col; i < first_col+cols_per_page; i++) 
          fprintf(outfile,"              %5d",i);

      fprintf(outfile,"\n               ");
      for(i=first_col; i < first_col+cols_per_page; i++) 
          fprintf(outfile,"          (%3d,%3d)",
                  Params->colorb[irrep][i][0], Params->colorb[irrep][i][1]);

      fprintf (outfile,"\n");
      for(i=0; i < rows; i++) {
          fprintf(outfile,"\n%5d  (%3d,%3d)",i,
                  Params->roworb[irrep][i][0], Params->roworb[irrep][i][1]);

          for(j=first_col; j < first_col+cols_per_page; j++)        
              fprintf (outfile,"%19.15f",matrix[i][j]);
        }

      fprintf (outfile,"\n");
    }

  /* Now print the remaining columns */
  if(last_page) {
      first_col = page*cols_per_page;

      fprintf(outfile,"\n           ");
      for(i=first_col; i < first_col+last_page; i++) 
	  fprintf(outfile,"              %5d",i);
      
      fprintf(outfile,"\n               ");
      for(i=first_col; i < first_col+last_page; i++) 
	  fprintf(outfile,"          (%3d,%3d)",
		  Params->colorb[irrep][i][0], Params->colorb[irrep][i][1]);

      fprintf (outfile,"\n");
      for(i=0; i < rows; i++) {
	  fprintf(outfile,"\n%5d  (%3d,%3d)",i,
		  Params->roworb[irrep][i][0], Params->roworb[irrep][i][1]);

	  for(j=first_col; j < first_col+last_page; j++)
	      fprintf (outfile,"%19.15f",matrix[i][j]);
	}

      fprintf (outfile,"\n");
    }

  return 0;

}
