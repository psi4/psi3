/*! \file 
    \ingroup (TRANSQT2)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

/* get_moinfo(): Routine to obtain basic orbital information from
** chkpt and compute the associated lookup arrays.
*/

namespace psi {
  namespace transqt2 {

    void get_moinfo(void)
    {
      int i, j, k, h, p, count;
      double escf;
      int *offset, this_offset;
      int *rstr_docc, *rstr_uocc, *tmparray, **ras_opi;
  
      chkpt_init(PSIO_OPEN_OLD);
      moinfo.nirreps = chkpt_rd_nirreps();
      moinfo.nmo = chkpt_rd_nmo();
      moinfo.nso = chkpt_rd_nso();
      moinfo.nao = chkpt_rd_nao();
      moinfo.labels = chkpt_rd_irr_labs();
      moinfo.enuc = chkpt_rd_enuc();
      escf = chkpt_rd_escf();
      moinfo.sopi = chkpt_rd_sopi();
      moinfo.mopi = chkpt_rd_orbspi();
      moinfo.clsdpi = chkpt_rd_clsdpi();
      moinfo.openpi = chkpt_rd_openpi();
      chkpt_close();

      moinfo.frdocc = get_frzcpi();
      moinfo.fruocc = get_frzvpi();

      moinfo.nfzc = moinfo.nfzv = 0;
      for(i=0; i < moinfo.nirreps; i++) {
	moinfo.nfzc += moinfo.frdocc[i];
	moinfo.nfzv += moinfo.fruocc[i];
      }

      /** Compute spatial-orbial reordering array(s) for the one-electron transformation **/
      if(ci_wfn(params.wfn)) {

	rstr_docc = init_int_array(moinfo.nirreps);
	rstr_uocc = init_int_array(moinfo.nirreps);
	ras_opi = init_int_matrix(4,moinfo.nirreps);

	// below, we need frdocc and fruocc as they appear in input, i.e., they
	// should not be zeroed out even if the frozen orbitals are to be 
	// transformed

	moinfo.pitz2corr_one = init_int_array(moinfo.nmo);
	if (!ras_set2(moinfo.nirreps, moinfo.nmo, 1, 1,
		      moinfo.mopi, moinfo.clsdpi, moinfo.openpi,
		      moinfo.frdocc, moinfo.fruocc,
		      rstr_docc, rstr_uocc,
		      ras_opi, moinfo.pitz2corr_one, 1, 0))
	  {
	    fprintf(outfile, "Error in ras_set().  Aborting.\n");
	    exit(1);
	  }
      }
      else {  /* CC only for now */
	if(params.ref == 0 || params.ref == 1) {
	  moinfo.pitz2corr_one = init_int_array(moinfo.nmo);
	  reorder_qt(moinfo.clsdpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc, 
		     moinfo.pitz2corr_one, moinfo.mopi, moinfo.nirreps);
	}
	else if(params.ref == 2) {
	  moinfo.pitz2corr_one_A = init_int_array(moinfo.nmo);
	  moinfo.pitz2corr_one_B = init_int_array(moinfo.nmo);
	  reorder_qt_uhf(moinfo.clsdpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc, 
			 moinfo.pitz2corr_one_A, moinfo.pitz2corr_one_B, moinfo.mopi, moinfo.nirreps);
	}
      }

      /* We want actpi to include only active orbitals */
      moinfo.actpi = init_int_array(moinfo.nirreps);
      for(h=0; h < moinfo.nirreps; h++)
	moinfo.actpi[h] = moinfo.mopi[h] - moinfo.frdocc[h] - moinfo.fruocc[h];

      /* SO and MO symmetry arrays */
      moinfo.sosym = init_int_array(moinfo.nso);
      for(h=0,count=0; h < moinfo.nirreps; h++)
	for(i=0; i < moinfo.sopi[h]; i++,count++)
	  moinfo.sosym[count] = h;

      moinfo.actsym = init_int_array(moinfo.nmo - moinfo.nfzc - moinfo.nfzv);
      for(h=0,count=0; h < moinfo.nirreps; h++)
	for(i=0; i < moinfo.actpi[h]; i++,count++)
	  moinfo.actsym[count] = h;

      moinfo.uoccpi = init_int_array(moinfo.nirreps);
      for(i=0; i < moinfo.nirreps; i++)
	moinfo.uoccpi[i] = moinfo.mopi[i] - moinfo.clsdpi[i] - moinfo.openpi[i];

      moinfo.nactive = moinfo.nmo - moinfo.nfzc - moinfo.nfzv;

      /** Compute spatial-orbial reordering array(s) for the two-electron transformation **/

      if(ci_wfn(params.wfn)) {
	// Now we need to translate the full Pitzer -> correlated mapping
	// array to one that involves only the active orbitals
	moinfo.pitz2corr_two = init_int_array(moinfo.nactive);
	for (h=0,j=0,k=0; h<moinfo.nirreps; h++) {
	  for (i=0; i<moinfo.mopi[h]; i++,j++) {
	    // j is now an absolute Pitzer MO index  
	    // k is an index for Pitzer orbitals NOT including frozen ones
	    if (i < moinfo.frdocc[h] || i >= (moinfo.mopi[h]-moinfo.fruocc[h])) 
	      continue;   
	    moinfo.pitz2corr_two[k++] = moinfo.pitz2corr_one[j] - moinfo.nfzc; 
	  }  
	}
	free_int_matrix(ras_opi);
	free(rstr_docc);
	free(rstr_uocc);
      }
      else { /* CC only */

	offset = init_int_array(moinfo.nirreps);
	for(h=1; h < moinfo.nirreps; h++)
	  offset[h] = offset[h-1] + moinfo.actpi[h-1];

	if(params.ref == 0 || params.ref == 1) {

	  moinfo.pitz2corr_two = init_int_array(moinfo.nactive);

	  count = 0;
	  for(h=0; h < moinfo.nirreps; h++) {
	    this_offset = offset[h];
	    for(p=0; p < moinfo.clsdpi[h] - moinfo.frdocc[h]; p++)
	      moinfo.pitz2corr_two[this_offset+p] = count++;
	  }
	  for(h=0; h < moinfo.nirreps; h++) {
	    this_offset = offset[h] + moinfo.clsdpi[h] - moinfo.frdocc[h];
	    for(p=0; p < moinfo.openpi[h]; p++)
	      moinfo.pitz2corr_two[this_offset+p] = count++;
	  }
	  for(h=0; h < moinfo.nirreps; h++) {
	    this_offset = offset[h] + moinfo.clsdpi[h] - moinfo.frdocc[h] + moinfo.openpi[h];
	    for(p=0; p < moinfo.uoccpi[h] - moinfo.fruocc[h]; p++)
	      moinfo.pitz2corr_two[this_offset+p] = count++;
	  }
	}
	else if(params.ref == 2) {
	  moinfo.pitz2corr_two_A = init_int_array(moinfo.nactive);
	  moinfo.pitz2corr_two_B = init_int_array(moinfo.nactive);

	  count = 0;
	  for(h=0; h < moinfo.nirreps; h++) {
	    this_offset = offset[h];
	    for(p=0; p < moinfo.clsdpi[h] - moinfo.frdocc[h] + moinfo.openpi[h]; p++)
	      moinfo.pitz2corr_two_A[this_offset+p] = count++;
	  }
	  for(h=0; h < moinfo.nirreps; h++) {
	    this_offset = offset[h] + moinfo.clsdpi[h] - moinfo.frdocc[h] + moinfo.openpi[h];
	    for(p=0; p < moinfo.uoccpi[h] - moinfo.fruocc[h]; p++)
	      moinfo.pitz2corr_two_A[this_offset+p] = count++;
	  }

	  count = 0;
	  for(h=0; h < moinfo.nirreps; h++) {
	    this_offset = offset[h];
	    for(p=0; p < moinfo.clsdpi[h] - moinfo.frdocc[h]; p++)
	      moinfo.pitz2corr_two_B[this_offset+p] = count++;
	  }
	  for(h=0; h < moinfo.nirreps; h++) {
	    this_offset = offset[h] + moinfo.clsdpi[h] - moinfo.frdocc[h];
	    for(p=0; p < moinfo.openpi[h] + moinfo.uoccpi[h] - moinfo.fruocc[h]; p++)
	      moinfo.pitz2corr_two_B[this_offset+p] = count++;
	  }
	}
	free(offset);

      }

      if(params.print_lvl) {
	fprintf(outfile,"\tChkpt Parameters:\n");
	fprintf(outfile,"\t--------------------\n");
	fprintf(outfile,"\tNumber of irreps     = %d\n",moinfo.nirreps);
	fprintf(outfile,"\tNumber of SOs        = %d\n",moinfo.nso);
	fprintf(outfile,"\tNumber of MOs        = %d\n",moinfo.nmo);
	fprintf(outfile,"\tNumber of active MOs = %d\n\n",moinfo.nactive);
	fprintf(outfile,
		"\tLabel\t# SOs\t# FZDC\t# DOCC\t# SOCC\t# VIRT\t# FZVR\n");
	fprintf(outfile,
		"\t-----\t-----\t------\t------\t------\t------\t------\n");
	for(i=0; i < moinfo.nirreps; i++) {
	  fprintf(outfile,
		  "\t %s\t   %d\t    %d\t    %d\t    %d\t    %d\t    %d\n",
		  moinfo.labels[i],moinfo.sopi[i],moinfo.frdocc[i],
		  moinfo.clsdpi[i]-moinfo.frdocc[i],moinfo.openpi[i],moinfo.uoccpi[i]-moinfo.fruocc[i],
		  moinfo.fruocc[i]);
	}
	fprintf(outfile,"\n\tNuclear Rep. energy (chkpt) =  %20.14f\n", moinfo.enuc);
	fprintf(outfile,  "\tSCF energy          (chkpt) =  %20.14f\n", escf);
      }
    }

    void cleanup(void)
    {
      int h;
      free(moinfo.sopi);
      free(moinfo.sosym);
      free(moinfo.mopi);
      free(moinfo.mosym);
      free(moinfo.actpi);
      free(moinfo.actsym);
      free(moinfo.clsdpi);
      free(moinfo.openpi);
      free(moinfo.uoccpi);
      free(moinfo.frdocc);
      free(moinfo.fruocc);
      for(h=0; h < moinfo.nirreps; h++)
	free(moinfo.labels[h]);
      free(moinfo.labels);
      if(params.ref == 0 || params.ref == 1) {
	free(moinfo.pitz2corr_one);
	free(moinfo.pitz2corr_two);
      }
      else if(params.ref == 2) {
	free(moinfo.pitz2corr_one_A);
	free(moinfo.pitz2corr_one_B);
	free(moinfo.pitz2corr_two_A);
	free(moinfo.pitz2corr_two_B);
      }
    }

  } // namespace transqt2
} // namespace psi
