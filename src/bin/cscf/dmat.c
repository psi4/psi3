/* $Log$
 * Revision 1.1  2000/02/04 22:52:29  evaleev
 * Initial revision
 *
/* Revision 1.6  1999/11/17 19:40:45  evaleev
/* Made all the adjustments necessary to have direct UHF working. Still doesn't work though..
/*
/* Revision 1.5  1999/11/02 23:23:42  evaleev
/* Commented out a line in dmat..
/*
/* Revision 1.4  1999/11/02 18:10:13  evaleev
/* Direct SCF improved
/*
/* Revision 1.3  1999/10/22 19:47:18  evaleev
/* A direct SCF-enabled version (set DIRECT_SCF=TRUE in input.dat).
/*
/* Revision 1.2  1999/08/17 19:04:14  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:25  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
/*
 * Revision 1.1  1991/06/15  20:22:21  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include <psio.h>
#include "includes.h"
#include "common.h"

void dmat()
{
   int i,j,k,l,ij,jj,kk,nn;
   int max, off, ntri;
   int ndocc,nsocc,nhocc;
   double ptempc,ptempo,ctmp;
   double *dmat, *dmato;
   extern double delta;
   struct symm *s;

   for (l=0; l < num_ir ; l++) {
      s = &scf_info[l];
      if (nn=s->num_so) {
         ndocc = s->nclosed;
         nsocc = s->nopen;
         nhocc = 0;
         if(s->nhalf) nhocc=1;
            
         for (i=ij=0; i < nn ; i++ ) {
            for (j=0; j < i; j++,ij++) {
               ptempc=ptempo=0.0;
               for (k=0; k < ndocc ; k++)
                  ptempc += 4.0*s->cmat[i][k]*s->cmat[j][k];

               for (k=ndocc; k < ndocc+nsocc ; k++)
                  ptempo += 2.0*s->occ_num[k]*s->cmat[i][k]*s->cmat[j][k];

               for (k=ndocc+nsocc; k < ndocc+nsocc+nhocc ; k++)
                  ptempo += 2.0*s->occ_num[k]*s->cmat[i][k]*s->cmat[j][k];

               if(iopen) {
                  s->dpmato[ij] = ptempo - s->pmato[ij];
                  s->pmato[ij] = ptempo;
                  }
               s->dpmat[ij] = ptempc+ptempo - s->pmat[ij];
               s->pmat[ij] = ptempc+ptempo;
               }
            ptempc=ptempo=0.0;
            for (k=0; k < ndocc ; k++) {
               ctmp=s->cmat[i][k];
               ptempc += 2.0*ctmp*ctmp;
               }
            for (k=ndocc; k < ndocc+nsocc ; k++) {
               ctmp=s->cmat[i][k];
               ptempo += ctmp*ctmp*s->occ_num[k];
               }
            for (k=ndocc+nsocc; k < ndocc+nsocc+nhocc ; k++) {
               ctmp=s->cmat[i][k];
               ptempo += ctmp*ctmp*s->occ_num[k];
               }
            if(iopen) {
               s->dpmato[ij] = ptempo - s->pmato[ij];
               s->pmato[ij] = ptempo;
               }
            s->dpmat[ij] = ptempc+ptempo - s->pmat[ij];
            s->pmat[ij] = ptempc+ptempo;
            ij++;
            }

         if(print & 4) {
            fprintf(outfile,
                       "\ntotal density matrix for irrep %s",s->irrep_label);
            print_array(s->pmat,nn,outfile);
            print_array(s->dpmat,nn,outfile);
            if(iopen) {
               fprintf(outfile,"\nopen-shell density matrix for irrep %s",
                                                              s->irrep_label);
               print_array(s->pmato,nn,outfile);
               print_array(s->dpmato,nn,outfile);
               }
            }
         }
      }

   /*------------------------
     Get full dpmat and dpmato
    ------------------------*/
   ntri = nbasis*(nbasis+1)/2;
   dmat = init_array(ntri);
   if (iopen)
     dmato = init_array(ntri);

   for(i=0;i<num_ir;i++) {
     max = scf_info[i].num_so;
     off = scf_info[i].ideg;
     for(j=0;j<max;j++) {
       jj = j + off;
       for(k=0;k<=j;k++) {
	 kk = k + off;
	 dmat[ioff[jj]+kk] = scf_info[i].dpmat[ioff[j]+k];
       }
     }
     if (iopen)
       for(j=0;j<max;j++) {
	 jj = j + off;
	 for(k=0;k<=j;k++) {
	   kk = k + off;
	   dmato[ioff[jj]+kk] = scf_info[i].dpmato[ioff[j]+k];
	 }
       }
   }

   if (direct_scf) {
     psio_open(itapDSCF, PSIO_OPEN_NEW);
     ctmp = 1.0;
     psio_write_entry(itapDSCF, "HF exchange contribution", (char *) &ctmp, sizeof(double));
     /*--- Decide what accuracy to request ---*/
     if (iconv <= 3)
       eri_cutoff = 1.0E-8;
     else if (iconv <= 5)
       eri_cutoff = 1.0E-11;
     else
       eri_cutoff = 1.0E-14;
     psio_write_entry(itapDSCF, "Integrals cutoff", (char *) &eri_cutoff, sizeof(double));
     psio_write_entry(itapDSCF, "Total SO Density", (char *) dmat, sizeof(double)*ntri);
     if (iopen) {
       psio_write_entry(itapDSCF, "Open-shell SO Density", (char *) dmato, sizeof(double)*ntri);
     }
     psio_close(itapDSCF, 1);
   }

   free(dmat);
   if (iopen)
     free(dmato);
   
}
