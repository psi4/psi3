#include <stdio.h>
#include <stdlib.h>
#include <psio.h>
#include <libciomr.h>
#include <iwl.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

void sort_oei(void)
{
  int h;
  int p,q,pq;
  int pnew, qnew;
  int psym, qsym;
  int nirreps;
  int *occ, *vir;
  int *cc_occ, *cc_vir;
  int *occ_sym, *vir_sym;
  int *occpi, *virtpi;
  int *occ_off, *vir_off;
  int docc;
  double *oei;
  psio_address next_hij, next_hab, next_hia;

  nirreps = moinfo.nirreps;
  occ = moinfo.occ; vir = moinfo.vir;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  cc_occ = moinfo.cc_occ; cc_vir = moinfo.cc_vir;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;

  /* Grab the frozen-core opertor */
  oei = init_array(moinfo.nactive*(moinfo.nactive+1)/2);
  iwl_rdone(params.FZCFile, oei, &moinfo.efzc, ioff, moinfo.nmo,
	    moinfo.nfzc, moinfo.nfzv, 0, 0, outfile);

  if(params.print_lvl > 5) {
      fprintf(outfile, "\n\tFrozen-Core Operator:\n");
      fprintf(outfile,   "\t---------------------");
      print_array(oei, moinfo.nactive, outfile);
    }

  hoo = (double ***) malloc(nirreps * sizeof(double **));
  hvv = (double ***) malloc(nirreps * sizeof(double **));
  hov = (double ***) malloc(nirreps * sizeof(double **));
  for(h=0; h < nirreps; h++) {
      hoo[h] = block_matrix(occpi[h],occpi[h]);
      hvv[h] = block_matrix(virtpi[h],virtpi[h]);
      hov[h] = block_matrix(occpi[h],virtpi[h]);
    }

  /* Filter out the frozen orbitals to generate the DPD-blocked
     one-electron integral lists */
  for(p=0; p < moinfo.nactive; p++) {
      for(q=0; q < moinfo.nactive; q++) {
	  pq = INDEX(p, q);

	  /* Check occ-occ class */
	  if(occ[p] && occ[q]) {
	      /* Get relative indices */
	      pnew = cc_occ[p];  qnew = cc_occ[q];
	      /* Get orbital symmetries */
	      psym = occ_sym[pnew]; qsym = occ_sym[qnew];
	      /* Shift symmetry-relative indices */
	      pnew -= occ_off[psym]; qnew -= occ_off[qsym];
	      /* Check orbital symmetry and put integral in place */
	      if(psym == qsym) hoo[psym][pnew][qnew] = oei[pq];
	    }

	  /* Check vir-vir class */
	  if(vir[p] && vir[q]) {
	      /* Get relative indices */
	      pnew = cc_vir[p]; qnew = cc_vir[q];
	      /* Get orbital symmetries */
	      psym = vir_sym[pnew]; qsym = vir_sym[qnew];
	      /* Shift symmetry-relative indices */
	      pnew -= vir_off[psym]; qnew -= vir_off[qsym];
	      /* Check orbital symmetry and put integral in place */
	      if(psym == qsym) hvv[psym][pnew][qnew] = oei[pq];
	    }

	  /* Check occ-vir class */
	  if(occ[p] && vir[q]) {
	      /* Get relative indices */
	      pnew = cc_occ[p]; qnew = cc_vir[q];
	      /* Get orbital symmetries */
	      psym = occ_sym[pnew]; qsym = vir_sym[qnew];
	      /* Shift symmetry-relative indices */
	      pnew -= occ_off[psym]; qnew -= vir_off[qsym];
	      /* Check orbital symmetry and put integral in place */
	      if(psym == qsym) hov[psym][pnew][qnew] = oei[pq];
	    }
	}
    }

  if(params.print_lvl > 2) {
      fprintf(outfile, "\n\thoo Matrix:\n");
      fprintf(outfile,   "\t-----------\n");
      for(h=0; h < nirreps; h++)
          print_mat(hoo[h], occpi[h], occpi[h], outfile);

      fprintf(outfile, "\n\thvv Matrix:\n");
      fprintf(outfile,   "\t-----------\n");
      for(h=0; h < nirreps; h++)
          print_mat(hvv[h], virtpi[h], virtpi[h], outfile);

      fprintf(outfile, "\n\thov Matrix:\n");
      fprintf(outfile,   "\t-----------\n");
      for(h=0; h < nirreps; h++)
          print_mat(hov[h], occpi[h], virtpi[h], outfile);
    }

  /* Dumping one-electron integrals in DPD format for use in ccdensity */
  next_hij = PSIO_ZERO;
  for(h=0; h < nirreps; h++) 
      if(occpi[h])
	  psio_write(CC_OEI, "h(i,j)", (char *) hoo[h][0],
		     occpi[h]*occpi[h]*sizeof(double), next_hij, &next_hij);

  next_hab = PSIO_ZERO;
  for(h=0; h < nirreps; h++) 
      if(virtpi[h])
	  psio_write(CC_OEI, "h(a,b)", (char *) hvv[h][0],
		     virtpi[h]*virtpi[h]*sizeof(double), next_hab, &next_hab);

  next_hia = PSIO_ZERO;
  for(h=0; h < nirreps; h++) 
      if(occpi[h] && virtpi[h])
	  psio_write(CC_OEI, "h(i,a)", (char *) hov[h][0],
		     occpi[h]*virtpi[h]*sizeof(double), next_hia,
		     &next_hia);

  free(oei);
}
