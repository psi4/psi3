#include <stdio.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#define EXTERN
#include "globals.h"

/*
** get_moinfo():  Routine to obtain basic orbital information from
** CHKPT and CC_INFO.
**
** T. Daniel Crawford, October 1996.
** Modified by TDC, March 1999.
*/

void get_moinfo(void)
{
  int i, j, h, errcod;
  int nactive;

  chkpt_init();
  moinfo.nirreps = chkpt_rd_nirreps();
  moinfo.nmo = chkpt_rd_nmo();
  moinfo.iopen = chkpt_rd_iopen();
  moinfo.labels = chkpt_rd_irr_labs();
  moinfo.enuc = chkpt_rd_enuc();
  moinfo.escf = chkpt_rd_escf();
  moinfo.orbspi = chkpt_rd_orbspi();
  moinfo.clsdpi = chkpt_rd_clsdpi();
  moinfo.openpi = chkpt_rd_openpi();
  chkpt_close();

  /* Get frozen and active orbital lookups from CC_INFO */
  moinfo.frdocc = init_int_array(moinfo.nirreps);
  moinfo.fruocc = init_int_array(moinfo.nirreps);
  psio_read_entry(CC_INFO, "Frozen Core Orbs Per Irrep",
		  (char *) moinfo.frdocc, sizeof(int)*moinfo.nirreps);
  psio_read_entry(CC_INFO, "Frozen Virt Orbs Per Irrep",
		  (char *) moinfo.fruocc, sizeof(int)*moinfo.nirreps);
  
  moinfo.occpi = init_int_array(moinfo.nirreps);
  moinfo.virtpi = init_int_array(moinfo.nirreps);
  psio_read_entry(CC_INFO, "Active Occ Orbs Per Irrep",
		  (char *) moinfo.occpi, sizeof(int)*moinfo.nirreps);
  psio_read_entry(CC_INFO, "Active Virt Orbs Per Irrep",
		  (char *) moinfo.virtpi, sizeof(int)*moinfo.nirreps);

  psio_read_entry(CC_INFO, "No. of Active Orbitals", (char *) &(nactive),
		  sizeof(int));
  moinfo.nactive = nactive;

  moinfo.occ_sym = init_int_array(nactive);
  moinfo.vir_sym = init_int_array(nactive);
  psio_read_entry(CC_INFO, "Active Occ Orb Symmetry",
		  (char *) moinfo.occ_sym, sizeof(int)*nactive);
  psio_read_entry(CC_INFO, "Active Virt Orb Symmetry",
		  (char *) moinfo.vir_sym, sizeof(int)*nactive);

  moinfo.occ_off = init_int_array(moinfo.nirreps);
  moinfo.vir_off = init_int_array(moinfo.nirreps);
  psio_read_entry(CC_INFO, "Active Occ Orb Offsets",
		  (char *) moinfo.occ_off, sizeof(int)*moinfo.nirreps);
  psio_read_entry(CC_INFO, "Active Virt Orb Offsets",
		  (char *) moinfo.vir_off, sizeof(int)*moinfo.nirreps);

  /* Adjust clsdpi array for frozen orbitals */
  for(i=0; i < moinfo.nirreps; i++)
      moinfo.clsdpi[i] -= moinfo.frdocc[i];

  moinfo.uoccpi = init_int_array(moinfo.nirreps);
  for(i=0; i < moinfo.nirreps; i++)
      moinfo.uoccpi[i] = moinfo.orbspi[i] - moinfo.clsdpi[i] -
			 moinfo.openpi[i] - moinfo.fruocc[i] -
			 moinfo.frdocc[i];

  moinfo.nfzc = moinfo.nfzv = moinfo.nclsd = moinfo.nopen = moinfo.nuocc = 0;
  for(h=0; h < moinfo.nirreps; h++) {
      moinfo.nfzc += moinfo.frdocc[h];
      moinfo.nfzv += moinfo.fruocc[h];
      moinfo.nclsd += moinfo.clsdpi[h];
      moinfo.nopen += moinfo.openpi[h];
      moinfo.nuocc += moinfo.uoccpi[h];
    }

  /* Get CC->QT and QT->CC active occupied and virtual reordering arrays */
  moinfo.qt_occ = init_int_array(nactive);
  moinfo.qt_vir = init_int_array(nactive);
  moinfo.cc_occ = init_int_array(nactive);
  moinfo.cc_vir = init_int_array(nactive);
  psio_read_entry(CC_INFO, "CC->QT Active Occ Order",
		   (char *) moinfo.qt_occ, sizeof(int)*nactive);
  psio_read_entry(CC_INFO, "CC->QT Active Virt Order",
		   (char *) moinfo.qt_vir, sizeof(int)*nactive);
  psio_read_entry(CC_INFO, "QT->CC Active Occ Order",
		   (char *) moinfo.cc_occ, sizeof(int)*nactive);
  psio_read_entry(CC_INFO, "QT->CC Active Virt Order",
		   (char *) moinfo.cc_vir, sizeof(int)*nactive);

  psio_read_entry(CC_INFO, "Reference Energy", (char *) &(moinfo.eref),
		  sizeof(double));
  psio_read_entry(CC_INFO, "CCSD Energy", (char *) &(moinfo.ecc),
		  sizeof(double));

  fprintf(outfile,"\n\tNuclear Rep. energy (chkpt)   = %20.15f\n",moinfo.enuc);
  fprintf(outfile,  "\tSCF energy          (chkpt)   = %20.15f\n",moinfo.escf);
  fprintf(outfile,  "\tReference energy    (file100) = %20.15f\n",moinfo.eref);
  fprintf(outfile,  "\tCCSD energy         (file100) = %20.15f\n",moinfo.ecc);
  fprintf(outfile,  "\tTotal CCSD energy   (file100) = %20.15f\n", 
          moinfo.eref+moinfo.ecc);
}

/* Frees memory allocated in get_moinfo(). */
void cleanup(void)
{
  int i;

  free(moinfo.orbspi);
  free(moinfo.clsdpi);
  free(moinfo.openpi);
  free(moinfo.uoccpi);
  free(moinfo.fruocc);
  free(moinfo.frdocc);
  for(i=0; i < moinfo.nirreps; i++)
      free(moinfo.labels[i]);
  free(moinfo.labels);
  free(moinfo.occ_sym);
  free(moinfo.vir_sym);
  free(moinfo.occ_off);
  free(moinfo.vir_off);
  free(moinfo.occpi);
  free(moinfo.virtpi);
  free(moinfo.qt_occ);
  free(moinfo.qt_vir);
}

