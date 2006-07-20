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

void get_moinfo(void)
{
  int i, h, p, count;
  int nfzc, nuocc, nopen, nclsd;
  double escf;
  int *offset, this_offset;
  
  chkpt_init(PSIO_OPEN_OLD);
  moinfo.nirreps = chkpt_rd_nirreps();
  moinfo.nmo = chkpt_rd_nmo();
  moinfo.nso = chkpt_rd_nso();
  moinfo.nao = chkpt_rd_nao();
  moinfo.labels = chkpt_rd_irr_labs();
  moinfo.enuc = chkpt_rd_enuc();
  escf = chkpt_rd_escf();
  moinfo.sopi = chkpt_rd_orbspi();
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

  /* MOs per irrep */
  moinfo.mopi = init_int_array(moinfo.nirreps);
  for(h=0; h < moinfo.nirreps; h++)
    moinfo.mopi[h] = moinfo.sopi[h] - moinfo.frdocc[h] - moinfo.fruocc[h];

  /* SO and MO symmetry arrays */
  moinfo.sosym = init_int_array(moinfo.nso);
  for(h=0,count=0; h < moinfo.nirreps; h++)
    for(i=0; i < moinfo.sopi[h]; i++,count++)
      moinfo.sosym[count] = h;

  moinfo.mosym = init_int_array(moinfo.nmo - moinfo.nfzc - moinfo.nfzv);
  for(h=0,count=0; h < moinfo.nirreps; h++)
    for(i=0; i < moinfo.mopi[h]; i++,count++)
      moinfo.mosym[count] = h;

  /* Compute spatial-orbial reordering arrays */
  if(params.ref == 0 || params.ref == 1) {
    moinfo.pitzer2qt = init_int_array(moinfo.nmo);
    moinfo.qt2pitzer = init_int_array(moinfo.nmo);
    reorder_qt(moinfo.clsdpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc, 
	       moinfo.pitzer2qt, moinfo.sopi, moinfo.nirreps);
    for(i=0; i < moinfo.nmo; i++)
      moinfo.qt2pitzer[moinfo.pitzer2qt[i]] = i;
  }
  else if(params.ref == 2) {
    moinfo.pitzer2qt_A = init_int_array(moinfo.nmo);
    moinfo.pitzer2qt_B = init_int_array(moinfo.nmo);
    moinfo.qt2pitzer_A = init_int_array(moinfo.nmo);
    moinfo.qt2pitzer_B = init_int_array(moinfo.nmo);
    reorder_qt_uhf(moinfo.clsdpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc, 
		   moinfo.pitzer2qt_A, moinfo.pitzer2qt_B, moinfo.sopi, moinfo.nirreps);
    for(i=0; i < moinfo.nmo; i++) {
      moinfo.qt2pitzer_A[moinfo.pitzer2qt_A[i]] = i;
      moinfo.qt2pitzer_B[moinfo.pitzer2qt_B[i]] = i;
    }
  }

  /* Adjust clsdpi array for frozen orbitals */
  for(i=0; i < moinfo.nirreps; i++)
    moinfo.clsdpi[i] -= moinfo.frdocc[i];

  moinfo.uoccpi = init_int_array(moinfo.nirreps);
  for(i=0; i < moinfo.nirreps; i++)
    moinfo.uoccpi[i] = moinfo.sopi[i] - moinfo.clsdpi[i] -
      moinfo.openpi[i] - moinfo.fruocc[i] -
      moinfo.frdocc[i];

  nclsd = nopen = nuocc = 0;
  for(i=0; i < moinfo.nirreps; i++) {
    nclsd += moinfo.clsdpi[i];
    nopen += moinfo.openpi[i];
    nuocc += moinfo.uoccpi[i];
  }
  nfzc = moinfo.nfzc;

  moinfo.nactive = nclsd + nopen + nuocc;

  /* more orbital reordering arrays */
  if(params.ref == 0 || params.ref == 1) {
    offset = init_int_array(moinfo.nirreps);
    for(h=1; h < moinfo.nirreps; h++)
      offset[h] = offset[h-1] + moinfo.mopi[h-1];

    moinfo.act2qt = init_int_array(moinfo.nactive);
    moinfo.qt2act = init_int_array(moinfo.nactive);

    count = 0;
    for(h=0; h < moinfo.nirreps; h++) {
      this_offset = offset[h];
      for(p=0; p < moinfo.clsdpi[h]; p++)
	moinfo.act2qt[this_offset+p] = count++;
    }
    for(h=0; h < moinfo.nirreps; h++) {
      this_offset = offset[h] + moinfo.clsdpi[h];
      for(p=0; p < moinfo.openpi[h]; p++)
	moinfo.act2qt[this_offset+p] = count++;
    }
    for(h=0; h < moinfo.nirreps; h++) {
      this_offset = offset[h] + moinfo.clsdpi[h] + moinfo.openpi[h];
      for(p=0; p < moinfo.uoccpi[h]; p++)
	moinfo.act2qt[this_offset+p] = count++;
    }
    free(offset);
  }
  else if(params.ref == 2) {

  }

  fprintf(outfile,"\n\tChkpt Parameters:\n");
  fprintf(outfile,"\t--------------------\n");
  fprintf(outfile,"\tNumber of irreps     = %d\n",moinfo.nirreps);
  fprintf(outfile,"\tNumber of MOs        = %d\n",moinfo.nmo);
  fprintf(outfile,"\tNumber of active MOs = %d\n\n",moinfo.nactive);
  fprintf(outfile,
	  "\tLabel\t# MOs\t# FZDC\t# DOCC\t# SOCC\t# VIRT\t# FZVR\n");
  fprintf(outfile,
	  "\t-----\t-----\t------\t------\t------\t------\t------\n");
  for(i=0; i < moinfo.nirreps; i++) {
    fprintf(outfile,
	    "\t %s\t   %d\t    %d\t    %d\t    %d\t    %d\t    %d\n",
	    moinfo.labels[i],moinfo.sopi[i],moinfo.frdocc[i],
	    moinfo.clsdpi[i],moinfo.openpi[i],moinfo.uoccpi[i],
	    moinfo.fruocc[i]);
  }
  fprintf(outfile,"\n\tNuclear Rep. energy (chkpt) =  %20.14f\n", moinfo.enuc);
  fprintf(outfile,  "\tSCF energy          (chkpt) =  %20.14f\n", escf);
}

void cleanup(void)
{
  int h;
  free(moinfo.sopi);
  free(moinfo.sosym);
  free(moinfo.mopi);
  free(moinfo.mosym);
  free(moinfo.clsdpi);
  free(moinfo.openpi);
  free(moinfo.uoccpi);
  free(moinfo.frdocc);
  free(moinfo.fruocc);
  for(h=0; h < moinfo.nirreps; h++)
    free(moinfo.labels[h]);
  free(moinfo.labels);
  if(params.ref == 0 || params.ref == 1) {
    free(moinfo.pitzer2qt);
    free(moinfo.qt2pitzer);
  }
  else if(params.ref == 2) {
    free(moinfo.pitzer2qt_A);
    free(moinfo.pitzer2qt_B);
    free(moinfo.qt2pitzer_A);
    free(moinfo.qt2pitzer_B);
  }
}
