#include <stdio.h>
#include <stdlib.h>
#include <psio.h>
#include <libciomr.h>
#include <iwl.h>
#include <dpd.h>
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
  dpdfile2 Hoo, Hov, Hvv;

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

  dpd_file2_init(&Hoo, CC_OEI, 0, 0, 0, "h(i,j)");
  dpd_file2_init(&Hvv, CC_OEI, 0, 1, 1, "h(a,b)");
  dpd_file2_init(&Hov, CC_OEI, 0, 0, 1, "h(i,a)");

  dpd_file2_mat_init(&Hoo);
  dpd_file2_mat_init(&Hvv);
  dpd_file2_mat_init(&Hov);

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
	      if(psym == qsym) {
		  Hoo.matrix[psym][pnew][qnew] = oei[pq];
		}
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
	      if(psym == qsym) {
		  Hvv.matrix[psym][pnew][qnew] = oei[pq];
		}
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
	      if(psym == qsym) {
		  Hov.matrix[psym][pnew][qnew] = oei[pq];
		}
	    }
	}
    }

  dpd_file2_mat_wrt(&Hoo);
  dpd_file2_mat_wrt(&Hvv);
  dpd_file2_mat_wrt(&Hov);

  dpd_file2_mat_close(&Hoo);
  dpd_file2_mat_close(&Hvv);
  dpd_file2_mat_close(&Hov);

  dpd_file2_close(&Hoo);
  dpd_file2_close(&Hvv);
  dpd_file2_close(&Hov);

  free(oei);
}
