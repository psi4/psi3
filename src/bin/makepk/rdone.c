
/* $Log$
 * Revision 1.1  2000/02/04 22:51:34  evaleev
 * Initial revision
 *
/* Revision 1.2  1997/08/25 21:52:31  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  22:06:29  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

void rdone(array)
   double *array;

   {
      int ilsti, nbuf;
      int ibufsz = 8942;
      int ibufs3 = 1491;
      int i, iijj, blk;
      int ior, ism, jor, jsm;


      union bufs {
         int *lbli;
         double *stvi;
         } buffer;

      buffer.stvi = (double *) init_array(ibufsz/2);

      do {
         sread(itap34,(char *) buffer.lbli,sizeof(int)*ibufsz);
         pos34 += sizeof(int)*ibufsz;
         pos34 = ((pos34-1+4096)/4096)*4096;
         ilsti=buffer.lbli[0];
         nbuf=buffer.lbli[1];

         if (print > 2) fprintf(outfile,"%5d\n",nbuf);

         for (i=0 ; i < nbuf ; i++) {
            jsm = buffer.lbli[i+2] >> 8;
            ior = jsm >> 3;
            ism = ior >> 8;
            ior = (ior & 255)-1;
            jsm = jsm & 7;
            jor = (buffer.lbli[i+2] & 255)-1;
            ior += ideg[ism];
            jor += ideg[jsm];
            iijj = ioff[ior]+jor;
            array[iijj]=buffer.stvi[i+ibufs3];
            }

          } while(!ilsti);
   free(buffer.stvi);
   }
