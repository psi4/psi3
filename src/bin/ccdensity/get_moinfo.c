#include <stdio.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <psio.h>
#include <file30.h>
#define EXTERN
#include "globals.h"

/*
** get_moinfo():  Routine to obtain basic orbital information from
** FILE30 and CC_INFO.
**
** T. Daniel Crawford, October 1996.
** Modified by TDC, March 1999.
*/

void get_moinfo(void)
{
  int i, j, h, errcod;
  int nactive;

  file30_init();
  moinfo.nirreps = file30_rd_nirreps();
  moinfo.nmo = file30_rd_nmo();
  moinfo.iopen = file30_rd_iopen();
  moinfo.labels = file30_rd_irr_labs();
  moinfo.enuc = file30_rd_enuc();
  moinfo.escf = file30_rd_escf();
  moinfo.orbspi = file30_rd_orbspi();
  moinfo.clsdpi = file30_rd_clsdpi();
  moinfo.openpi = file30_rd_openpi();
  file30_close();

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

  /* Get CC->QT active occupied and virtual reordering arrays */
  moinfo.qt_occ = init_int_array(nactive);
  moinfo.qt_vir = init_int_array(nactive);
  psio_read_entry(CC_INFO, "CC->QT Active Occ Order",
		   (char *) moinfo.qt_occ, sizeof(int)*nactive);
  psio_read_entry(CC_INFO, "CC->QT Active Virt Order",
		   (char *) moinfo.qt_vir, sizeof(int)*nactive);

  psio_read_entry(CC_INFO, "Reference Energy", (char *) &(moinfo.eref),
		  sizeof(double));
  psio_read_entry(CC_INFO, "CCSD Energy", (char *) &(moinfo.ecc),
		  sizeof(double));

  fprintf(outfile,"\n\tNuclear Rep. energy (file30)  = %20.15f\n",moinfo.enuc);
  fprintf(outfile,  "\tSCF energy          (file30)  = %20.15f\n",moinfo.escf);
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

