
/* $Log$
 * Revision 1.1  2000/02/04 22:50:48  evaleev
 * Initial revision
 *
/* Revision 1.2  1997/08/25 21:53:44  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  22:45:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

famat_c()
   {

      int i,j;
      double *hh,*tt,*ff;

      hh = (double *) init_array(nbstri);
      tt = (double *) init_array(nbstri);
      ff = (double *) init_array(nbstri);

      srew(work);

      for(i=0; i < natom3 ; i++) {
         rread(itap44,(char *) hh,sizeof(double)*nbstri,ha_loc[i]);
         sread(work,(char *) tt,sizeof(double)*nbstri);
         for(j=0; j < nbstri ; j++)
            ff[j] = hh[j]+tt[j];
         rwrit(itap44,(char *) ff,sizeof(double)*nbstri,fa_loc[i]);

         if(print & 8) {
            fprintf(outfile,"\nfa matrix i = %5d\n",i);
            print_array(ff,nbfso,outfile);
            }
         }

      free(hh);
      free(tt);
      free(ff);
      fflush(outfile);
      }
