/* $Log$
 * Revision 1.1  2000/02/04 22:52:31  evaleev
 * Initial revision
 *
/* Revision 1.6  1999/11/11 16:00:38  localpsi
/* Fixed minor bug in occupations.  STB (11/11/99)
/*
/* Revision 1.5  1999/11/04 19:24:29  localpsi
/* STB (11/4/99) - Added the orb_mix feature which is equivalent to guess = mix
/* in G94 and also fixed restarting so that if you have different wavefuntions,
/* everything works.  Also if you specify no DOCC and SOCC and restart, if the
/* wavefunctions are different, it will guess again.
/*
/* Revision 1.4  1999/11/02 18:10:13  evaleev
/* Direct SCF improved
/*
/* Revision 1.3  1999/08/17 19:04:15  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.2  1999/08/11 19:24:53  evaleev
/* Unhardwired the size of the ioff array (set it to 1024 for now) and increased MAX_BASIS to 1024.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:27  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"
#include <file30.h>

init_scf()
{
   int i,jj;
   int nn,isadr;
   int nkind,junk;
   PSI_FPTR next;
   int degen[20],*num_so;
   double *real_dum;
   char char_dum[80];
   char **irr_labs;

   i10 = (int *) init_int_array(200);
   real_dum = (double *) init_array(MAX_BASIS);

   wreadw(itap30,(char *) i10,sizeof(int)*200,(PSI_FPTR) sizeof(int)*100,&next);

   ioff[0] = 0;
   for (i = 1; i < MAXIOFF ; i++) {
      ioff[i] = ioff[i-1] + i;
      }

/* read header information from integral tape */
/* EFV 10/24/95 I left this junk here for now */

   if (!use_iwl) {
     num_so = init_int_array(20);
     wreadw(itap34,(char *) (&nkind),sizeof(int)*1,(PSI_FPTR) 0,&next);
     wreadw(itap34,(char *) (&junk),sizeof(int)*1,next,&next);
     wreadw(itap34,(char *) (char_dum),sizeof(char)*80,next,&next);
     wreadw(itap34,(char *) (&repnuc),sizeof(double)*1,next,&next);
     wreadw(itap34,(char *) (&num_ir),sizeof(int)*1,next,&next);
     wreadw(itap34,(char *) (degen),sizeof(int)*num_ir,next,&next);
     wreadw(itap34,(char *) (char_dum),sizeof(int)*num_ir,next,&next);
     wreadw(itap34,(char *) (num_so),sizeof(int)*num_ir,next,&next);
     wreadw(itap34,(char *) (&junk),sizeof(int)*1,next,&next);
     wreadw(itap34,(char *) (real_dum),sizeof(int)*2*junk,next,&next);
     wreadw(itap34,(char *) (&junk),sizeof(int)*1,next,&next);
     wreadw(itap34,(char *) (real_dum),sizeof(int)*junk,next,&next);
     wreadw(itap34,(char *) (real_dum),sizeof(int)*junk,next,&next);
     /* set integral file pointer to sector boundary */
     isadr = i2sec(next) + 1;
     rsetsa(itap34,isadr);
     pos34 = sec2i(--isadr);
     free(num_so);

     if (nkind != 1 && nkind != 2) {
        fprintf(outfile,"integral file screwed up, fix ints somebody!!!\n");
        /*exit(1);*/
     }
   }

/* EFV 10/24/98 All requests for file30 should be handled with libfile30
   but for now I'll use wreadw */
   num_ir = file30_rd_nirreps();
   num_so = file30_rd_sopi();
   repnuc = file30_rd_enuc();
   irr_labs = file30_rd_irr_labs();


/* now initialize scf_info */
   
   n_so_typs=0;
   mxcoef=0;
   mxcoef2=0;
   nsfmax=0;
   nbasis=0;
   
   scf_info = (struct symm *) malloc(sizeof(struct symm)*num_ir);

   jj=0;
   for(i=0; i < num_ir ; i++) {
      scf_info[i].num_so = nn = num_so[i];
/* EFV 10/24/98 degeneracy = 1
      scf_info[i].degeneracy = degen[i]; */
      scf_info[i].nclosed = 0;
      scf_info[i].nopen = 0;
      scf_info[i].nhalf = 0;
      scf_info[i].os_num = 0;
      scf_info[i].ideg = 0;

      scf_info[i].irrep_label = irr_labs[i];
/*      scf_info[i].irrep_label[4] = '\0';*/
      jj += 4;

      nbasis += nn;
      if (nn) {
         n_so_typs++;
         if (nn > nsfmax) nsfmax = nn;
         mxcoef += nn*nn;
         mxcoef2 += ioff[nn];

         scf_info[i].smat = (double *) init_array(ioff[nn]);
         scf_info[i].tmat = (double *) init_array(ioff[nn]);
         scf_info[i].hmat = (double *) init_array(ioff[nn]);
         scf_info[i].pmat = (double *) init_array(ioff[nn]);
         scf_info[i].pmato = (double *) NULL;
         scf_info[i].pmat2 = (double *) NULL;
         scf_info[i].pmato2 = (double *) NULL;
         scf_info[i].dpmat = (double *) init_array(ioff[nn]);
         scf_info[i].dpmato = (double *) NULL;
         scf_info[i].fock_pac = (double *) init_array(ioff[nn]);
         scf_info[i].fock_eff = (double *) NULL;
         scf_info[i].fock_open = (double *) NULL;
         scf_info[i].gmat = (double *) init_array(ioff[nn]);
         scf_info[i].gmato = (double *) NULL;
         scf_info[i].occ_num = (double *) init_array(nn);
         scf_info[i].fock_evals = (double *) init_array(nn);
         scf_info[i].cmat = (double **) init_matrix(nn,nn);
	 /* TDC(6/19/96) - Added array for saving original MO vector */
	 scf_info[i].cmat_orig = (double **) init_matrix(nn,nn);
         scf_info[i].sahalf = (double **) init_matrix(nn,nn);
         /* STB(4/1/98) - Added array to store the eigenvalues of the
	                  core hamiltonian for mo guessing*/
         scf_info[i].hevals = (double *) init_array(nn);
         }
     }
   /* read in number of atoms and nuclear charges and total number of MO*/
   natom = file30_rd_natom();
   zvals = file30_rd_zvals();
   nbfso = file30_rd_nso();
   
/* Initialize arrays to hold energy and symmetry arrays */
   ener_tot = (double *) init_array(nbfso);
   symm_tot = (int *) init_int_array(nbfso);
   
} 

init_scf2()
{
   int i,j,k;
   int n,nn,m,mm;
   int junk;
   int opconst,outbuf,mxcoef3,ntri,mtri;
   struct symm *s;

   opconst = (iopen) ? 3 : 2;

   mtri = ioff[scf_info[0].num_so];
   mxcoef3 = opconst*(mtri*(mtri+1)/2);

   scf_info[0].ideg = 0;

   so2symblk = init_int_array(nbasis);

   junk=j=0;
   for (i=0; i < num_ir ; i++) {
      s = &scf_info[i];
      nn = s->num_so;
      if (i) {
         m = ioff[nn];
         mm = m*(m+1)/2;
         mxcoef3 += opconst*ioff[nn]*mtri;
         mxcoef3 += opconst*mm;
         mtri += m;
         if (nn <= 0) s->ideg = scf_info[i-1].ideg;
         else {
            do {
               n=scf_info[j].num_so;
               j++;
               } while(!n); 
            s->ideg = scf_info[i-1].ideg+n;
            }
         }
       for(k=s->ideg;k<s->ideg+nn;k++)
         so2symblk[k] = i;
      if(s->nopen || s->nhalf) s->os_num = junk++;
      }
   ntri = nbasis*(nbasis+1)/2;

   if (direct_scf) { /* No I/O will be done */
     fprintf(outfile,
	     "\n  direct formation of the Fock matrix is requested\n");
   }
   else {            /* Figure out the size of the PK-buffers */
     readflg = 0;
     maxbuf = mxcoef3/opconst+5;
     outbuf = 884736;
     outbuf = MIN0(outbuf,mxcoef3);
     if(outbuf == 884736) {
       int pass = mxcoef3/outbuf+1;
       readflg = 1;
       maxbuf = (mxcoef3/pass)/opconst + 2;
       if(iopen) maxbuf /= 2;
#if defined(AIXV3)||defined(SGI)
       maxbuf= (iopen) ? 8192*3 : 8192*5;
       pass = mxcoef3/(maxbuf*opconst)+1;
#endif
       fprintf(outfile,
	       "\n  using buffered io, %d buffers, each %d bytes in size\n",
	       pass,maxbuf*opconst*sizeof(double));
       if(print) fprintf(outfile,"  mxcoef3 = %d maxbuf = %d\n",mxcoef3,maxbuf);
       if(print) fprintf(outfile,"  outbuf = %d\n",outbuf);
     }
     else fprintf(outfile,"\n  keeping integrals in %d bytes of core\n",
		  maxbuf*opconst*sizeof(double));
   }
     
   fflush(outfile);
   return;
   
}
