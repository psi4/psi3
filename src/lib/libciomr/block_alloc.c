
/* $Log$
 * Revision 1.2  2003/05/22 06:20:06  crawdad
 * Corrected most of the libraries and modules to use proper PSI_RETURN_XX
 * values from psifiles.h.  Modified ccdensity, ccenergy, cchbar, cclambda,
 * ccsort, cctriples, cis, cphf, cusp, localize, stable, libchkpt, libciomr,
 * libdpd, libipv1, libpsio, libqt, and tocprint.
 * -TDC
 *
/* Revision 1.1.1.1  2000/02/04 22:53:17  evaleev
/* Started PSI 3 repository
/*
/* Revision 2.3  1994/06/02 02:30:26  seidl
/* add test for too many rows
/*
 * Revision 1.1.1.1  1994/05/02  17:04:26  cljanss
 * The May 1, 1994 version of psi as on the CCQC machines.
 *
 * Revision 2.2  1991/06/15  18:53:37  seidl
 * remove include of common.h and definition of EXTERN
 *
 * Revision 2.1  1991/06/15  18:28:44  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id$";

#include <psifiles.h>
#include "includes.h"

double *** block_mat_alloc(n_so_typs,num_ir,num_so)
   int n_so_typs,num_ir,*num_so;

   {
      double ***array;
      int i,blk,nn;

      if ((array = (double ***) malloc(sizeof(double **)*n_so_typs))==NULL) {
          fprintf(stderr,"trouble in block_mat_alloc\n");
          exit(PSI_RETURN_FAILURE);
          }

      for (i=0,blk=0; i < num_ir ; i++) {
         if (nn=num_so[i]) {
            array[blk] = (double **) init_matrix(nn,nn);
            blk++;
            }
         }

      return(array);
      }

void block_mat_dealloc(array,num_ir,num_so)
   double ***array;
   int num_ir, *num_so;

   {
      int i;
      int blk=0;

      for (i=0; i < num_ir ; i++) {
         if (num_so[i]) {
            free_matrix(array[blk],num_so[i]);
            blk++;
            }
         }

      free(array);
      }

double ** block_arr_alloc(n_so_typs,num_ir,num_so)
   int n_so_typs,num_ir,*num_so;

   {
      double **array;
      int i;
      int j=0;

      if ((array = (double **) malloc(sizeof(double *)*n_so_typs))==NULL) {
          fprintf(stderr,"trouble in block_arr_alloc\n");
          exit(PSI_RETURN_FAILURE);
          }

      for (i=0; i < num_ir ; i++) {
         if (num_so[i]) {
            int nget = num_so[i]*(num_so[i]+1)/2;
            if (j>=n_so_typs) {
                fprintf(stderr,"block_arr_alloc: too many rows\n");
                exit(PSI_RETURN_FAILURE);
              }
            array[j] = (double *) init_array(nget);
            j++;
            }
         }
      return(array);
      }

void block_arr_dealloc(array,n_so_typs)
   double **array;

   {
      int i;

      for (i=0; i < n_so_typs ; i++) {
         free(array[i]);
         }

      free(array);
      }

