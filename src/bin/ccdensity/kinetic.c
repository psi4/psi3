#include <stdio.h>
#include <libciomr.h>
#include <iwl.h>
#include <file30.h>
#include <dpd.h>
#include <qt.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

#define IOFF_MAX 32641
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

void kinetic(void)
{
  int nmo, noei, stat, i, I, h, j, nclsd;
  int *order, *doccpi, *ioff;
  double junk, tcorr, vcorr, tref, vref, ttot, vtot;
  double *s, *t, **T, **S, **scf_pitzer, **scf_qt, **X;

  /*** Build ioff ***/
  ioff = init_int_array(IOFF_MAX);
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) ioff[i] = ioff[i-1] + i;

  nmo = moinfo.nmo;
  noei = nmo * (nmo + 1)/2;

  /*** Get the Pitzer -> QT reordering array ***/
  order = init_int_array(nmo);

  /* doccpi array must include frozen orbitals for reorder_qt() */
  doccpi = init_int_array(moinfo.nirreps);
  for(h=0; h < moinfo.nirreps; h++) 
      doccpi[h] = moinfo.frdocc[h] + moinfo.clsdpi[h];

  reorder_qt(doccpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc, 
             order, moinfo.orbspi, moinfo.nirreps);

  /*** Reorder the SCF eigenvectors to QT ordering */
  file30_init();
  scf_pitzer = file30_rd_scf();
  file30_close();

  scf_qt = block_matrix(nmo, nmo);
  for(i=0; i < nmo; i++) {
      I = order[i];  /* Pitzer --> QT */
      for(j=0; j < nmo; j++) scf_qt[j][I] = scf_pitzer[j][i];
    }

  /*
  fprintf(outfile, "\tSCF eigenvectors (Pitzer):\n");
  print_mat(scf_pitzer,nmo,nmo,outfile);
  fprintf(outfile, "\tSCF eigenvectors (QT):\n");
  print_mat(scf_qt,nmo,nmo,outfile);
  */

  /*** Transform the kinetic energy integrals to the MO basis ***/

  t = init_array(noei);
  stat = iwl_rdone(PSIF_OEI,PSIF_SO_T,t,noei,0,0,outfile);
  s = init_array(noei);
  stat = iwl_rdone(PSIF_OEI,PSIF_SO_S,s,noei,0,0,outfile);

  T = block_matrix(nmo,nmo);
  S = block_matrix(nmo,nmo);
  for(i=0; i < nmo; i++)
      for(j=0; j < nmo; j++) {
          T[i][j] = t[INDEX(i,j)];
          S[i][j] = s[INDEX(i,j)];
        }

  /*
  fprintf(outfile, "\tKinetic energy integrals (SO):\n");
  print_mat(T,nmo,nmo,outfile);
  fprintf(outfile, "\tOverlap integrals (SO):\n");
  print_mat(S,nmo,nmo,outfile);
  */

  X = block_matrix(nmo,nmo);

  C_DGEMM('t','n',nmo,nmo,nmo,1,&(scf_qt[0][0]),nmo,&(T[0][0]),nmo,
          0,&(X[0][0]),nmo);
  C_DGEMM('n','n',nmo,nmo,nmo,1,&(X[0][0]),nmo,&(scf_qt[0][0]),nmo,
          0,&(T[0][0]),nmo);

  /*
  fprintf(outfile, "\tKinetic energy integrals (MO):\n");
  print_mat(T,nmo,nmo,outfile);
  fprintf(outfile, "\tOne-particle density:\n");
  print_mat(moinfo.opdm,nmo,nmo,outfile);
  */

  /*** Contract the correlated kinetic energy ***/

  tcorr = 0.0;
  for(i=0; i < nmo; i++)
      for(j=0; j < nmo; j++)
          tcorr += T[i][j] * moinfo.opdm[i][j];

  /*** Compute the SCF kinetic energy ***/
  
  tref = 0.0;
  nclsd = moinfo.nfzc + moinfo.nclsd;
  for(i=0; i < nclsd; i++)
      tref += T[i][i] * 2;
  for(i=nclsd; i < nclsd+moinfo.nopen; i++)
      tref += T[i][i];

  /*** Compute the virial ratios ***/
  ttot = tcorr + tref;
  vtot = moinfo.eref + moinfo.ecc - ttot;
  vref = moinfo.eref - tref;
  vcorr = moinfo.ecc - tcorr;
  
  fprintf(outfile,"\n\tVirial Theorem Data:\n");
  fprintf(outfile,  "\t--------------------\n");
  fprintf(outfile,"\tKinetic energy (ref)   = %20.15f\n", tref);
  fprintf(outfile,"\tKinetic energy (corr)  = %20.15f\n", tcorr);
  fprintf(outfile,"\tKinetic energy (total) = %20.15f\n", ttot);

  fprintf(outfile,"\t-V/T (ref)             = %20.15f\n", -vref/tref);
  fprintf(outfile,"\t-V/T (corr)            = %20.15f\n", -vcorr/tcorr);
  fprintf(outfile,"\t-V/T (total)           = %20.15f\n", -vtot/ttot);

  fflush(outfile);

  /*** Release memory ***/
  free_block(X);
  free_block(T);
  free(t);
  free_block(scf_qt);
  free_block(scf_pitzer);
  free(doccpi);
  free(order);
  free(ioff);
}
