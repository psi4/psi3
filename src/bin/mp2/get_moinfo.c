#include <stdio.h>
#include <stdlib.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <ccfiles.h>
#include "moinfo.h"
#include "params.h"
#define EXTERN
#include "globals.h"

void get_moinfo(void)
{
  int i;
  
  chkpt_init(PSIO_OPEN_OLD);

  mo.nmo = chkpt_rd_nmo();
  mo.nso = chkpt_rd_nso();
  mo.nirreps = chkpt_rd_nirreps();
  mo.irreplabels = chkpt_rd_irr_labs();

  mo.mopi = chkpt_rd_orbspi();
  mo.doccpi = chkpt_rd_clsdpi();
  
  mo.Enuc = chkpt_rd_enuc();
  mo.Escf = chkpt_rd_escf();
  mo.scfevals = chkpt_rd_evals();
  
  mo.fzdoccpi = get_frzcpi();
  mo.fzvirtpi = get_frzvpi();
  
  chkpt_close();

  mo.nfzdocc = 0;
  mo.nfzvirt = 0;
  for (i=0; i<mo.nirreps; i++) {
    mo.nfzdocc += mo.fzdoccpi[i];
    mo.nfzvirt += mo.fzvirtpi[i];
  }

  mo.virtpi = init_int_array(mo.nirreps);
  for(i=0; i < mo.nirreps; i++) {
    mo.virtpi[i] = mo.mopi[i]-mo.doccpi[i];
  }
  
  mo.ndocc = 0;
  mo.nvirt = 0;
  for(i=0; i < mo.nirreps; i++) {
    mo.ndocc += mo.doccpi[i];
    mo.nvirt += mo.mopi[i] - mo.doccpi[i];
  }
  
  mo.actdoccpi = init_int_array(mo.nirreps);
  mo.actvirtpi = init_int_array(mo.nirreps);
  for(i=0; i < mo.nirreps; i++) {
    mo.actdoccpi[i] = mo.doccpi[i]-mo.fzdoccpi[i];
    mo.actvirtpi[i] = mo.virtpi[i]-mo.fzvirtpi[i];
  }

  mo.nactdocc = 0;
  mo.nactvirt = 0;
  for (i=0; i < mo.nirreps; i++) {
    mo.nactdocc += mo.actdoccpi[i];
    mo.nactvirt += mo.actvirtpi[i];
  }
 
  mo.nactmo = mo.nactdocc + mo.nactvirt;
  
  mo.actdoccsym = init_int_array(mo.nactmo);
  mo.actvirtsym = init_int_array(mo.nactmo);
  psio_read_entry(CC_INFO, "Active Occ Orb Symmetry",
	         (char *) mo.actdoccsym, sizeof(int)*mo.nactmo);
  psio_read_entry(CC_INFO, "Active Virt Orb Symmetry",
		 (char *) mo.actvirtsym, sizeof(int)*mo.nactmo);

  mo.docc_off = init_int_array(mo.nirreps);
  mo.virt_off = init_int_array(mo.nirreps);
  psio_read_entry(CC_INFO, "Active Occ Orb Offsets",
	          (char *) mo.docc_off, sizeof(int)*mo.nirreps);
  psio_read_entry(CC_INFO, "Active Virt Orb Offsets",
		  (char *) mo.virt_off, sizeof(int)*mo.nirreps);

  mo.qt_docc = init_int_array(mo.nactmo);
  mo.qt_virt = init_int_array(mo.nactmo);

  psio_read_entry(CC_INFO, "CC->QT Active Occ Order",
                 (char *) mo.qt_docc, sizeof(int)*mo.nactmo);
  psio_read_entry(CC_INFO, "CC->QT Active Virt Order",
	         (char *) mo.qt_virt, sizeof(int)*mo.nactmo);
	      
  fprintf(outfile,"\n");
  fprintf(outfile,"\tChkpt Parameters:\n");
  fprintf(outfile,"\t--------------------\n");
  fprintf(outfile,"\tNumber of irreps     = %d\n",mo.nirreps);
  fprintf(outfile,"\tNumber of MOs        = %d\n",mo.nmo);
  fprintf(outfile,"\n");
  fprintf(outfile,
    "\tLabel\tFZDC\tACTD\tDOCC\tACTV\tFZVI\tVIRT\tMOs\n");
  fprintf(outfile,
    "\t-----\t----\t----\t----\t----\t----\t----\t---\n");
  for(i=0; i < mo.nirreps; i++) {
    fprintf(outfile,
    "\t  %s \t  %d\t  %d\t  %d\t  %d\t  %d\t  %d\t %d\n",
	    mo.irreplabels[i],mo.fzdoccpi[i],mo.actdoccpi[i],mo.doccpi[i],
	    mo.actvirtpi[i],mo.fzvirtpi[i],mo.virtpi[i],mo.mopi[i]);
  }
  
  fprintf(outfile,"\n");
  fprintf(outfile,"\tNuclear Rep. energy \t=\t  %.12f\n",mo.Enuc);
  fprintf(outfile,"\tSCF energy          \t=\t%.12f\n",mo.Escf);
}
