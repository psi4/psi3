/* $Log$
 * Revision 1.1  2000/02/04 22:52:31  evaleev
 * Initial revision
 *
/* Revision 1.3  1999/11/02 18:10:14  evaleev
/* Direct SCF improved
/*
/* Revision 1.2  1999/08/17 19:04:16  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:27  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"
#include <iwl.h>

void rdone(oei)
   int oei;

{
   int ilsti, nbuf;
   int ibufsz = 8942;
   int ibufs3 = 1491;
   int i, iijj;
   int ior, ism, jor, jsm;

   union bufs {
      int *lbli;
      double *stvi;
      } buffer;

   buffer.stvi = (double *) init_array(ibufsz);

   do {
      sread(itap34,(char *) (buffer.lbli),sizeof(int)*ibufsz);
      pos34 += sizeof(int)*ibufsz;
      pos34 = ((pos34-1+4096)/4096)*4096;
      ilsti=buffer.lbli[0];
      nbuf=buffer.lbli[1];

      for (i=0 ; i < nbuf ; i++) {
         jsm = buffer.lbli[i+2] >> 8;
         ior = jsm >> 3;
         ism = ior >> 8;
         ior = (ior & 255)-1;
         jsm = jsm & 7;
         jor = (buffer.lbli[i+2] & 255)-1;
         iijj = ioff[ior]+jor;
         switch(oei) {
         case SMAT:
            scf_info[ism].smat[iijj]=buffer.stvi[i+ibufs3];
            break;
         case TMAT:
            scf_info[ism].tmat[iijj]=buffer.stvi[i+ibufs3];
            break;
         case VMAT:
            scf_info[ism].hmat[iijj]=
                               buffer.stvi[i+ibufs3]+scf_info[ism].tmat[iijj];
            break;
            }
         }
       } while(!ilsti);

   free(buffer.stvi);
   }


void rdone_iwl()
{
  int stat;
  int ntri = ioff[nbasis];
  int i,j,k,jj,kk;
  int max,off;
  double *ints;
  double e_fzc;

  /* If it's a direct SCF run - tell CINTS to compute one-electron integrals */
  if (direct_scf) {
    stat = system("cints --oeints");
    switch (stat) {
    case 0:
      /* CINTS ran successfully - continue */
      break;

    default:
      /* Something went wrong */
      fprintf(outfile,"  rdone_iwl: System call to CINTS failed. Check to see if it's in your PATH\n");
      fprintf(stderr,"System call to CINTS failed. Check to see if it's in your PATH.\n");
      exit(1);
    }
  }

  ints = init_array(ntri);

  /* S integrals */
  stat = iwl_rdone_all(itapS,ntri,ints,&e_fzc,delete_1e);
  for(i=0;i<num_ir;i++) {
    max = scf_info[i].num_so;
    off = scf_info[i].ideg;
    for(j=0;j<max;j++) {
      jj = j + off;
      for(k=0;k<=j;k++) {
	kk = k + off;
	scf_info[i].smat[ioff[j]+k] = ints[ioff[jj]+kk];
      }
    }
  }

  /* T integrals */
  stat = iwl_rdone_all(itapT,ntri,ints,&e_fzc,delete_1e);
  for(i=0;i<num_ir;i++) {
    max = scf_info[i].num_so;
    off = scf_info[i].ideg;
    for(j=0;j<max;j++) {
      jj = j + off;
      for(k=0;k<=j;k++) {
	kk = k + off;
	scf_info[i].tmat[ioff[j]+k] = ints[ioff[jj]+kk];
      }
    }
  }

  /* V integrals */
  stat = iwl_rdone_all(itapV,ntri,ints,&e_fzc,delete_1e);
  for(i=0;i<num_ir;i++) {
    max = scf_info[i].num_so;
    off = scf_info[i].ideg;
    for(j=0;j<max;j++) {
      jj = j + off;
      for(k=0;k<=j;k++) {
	kk = k + off;
	scf_info[i].hmat[ioff[j]+k] = ints[ioff[jj]+kk] + scf_info[i].tmat[ioff[j]+k];
      }
    }
  }
  free(ints);
}

