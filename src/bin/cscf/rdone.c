/* $Log$
 * Revision 1.2  2000/10/13 19:51:21  evaleev
 * Cleaned up a lot of stuff in order to get CSCF working with the new "Mo-projection-capable" INPUT.
 *
/* Revision 1.1.1.1  2000/02/04 22:52:31  evaleev
/* Started PSI 3 repository
/*
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

