#include <stdio.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <file30.h>
#include <psio.h>
#define EXTERN
#include "globals.h"

/*
** get_moinfo():  Routine to obtain basic orbital information from
** FILE30 and CC_INFO.
**
** T. Daniel Crawford, October 1996
** Modified by TDC, March 1999
*/

void get_moinfo(void)
{
  int i, h, p, q, errcod, nactive, nirreps;
  double ***C;
  psio_address next;

  file30_init();
  moinfo.nirreps = file30_rd_nirreps();
  moinfo.nmo = file30_rd_nmo();
  moinfo.nso = file30_rd_nso();
  moinfo.iopen = file30_rd_iopen();
  moinfo.labels = file30_rd_irr_labs();
  moinfo.enuc = file30_rd_enuc();
  moinfo.escf = file30_rd_escf();
  moinfo.orbspi = file30_rd_orbspi();
  moinfo.clsdpi = file30_rd_clsdpi();
  moinfo.openpi = file30_rd_openpi();
  moinfo.phase = file30_rd_phase_check();
  file30_close();

  nirreps = moinfo.nirreps;

  psio_read_entry(CC_INFO, "Reference Wavefunction", (char *) &(params.ref), 
		  sizeof(int));

  /* Get frozen and active orbital lookups from CC_INFO */
  moinfo.frdocc = init_int_array(nirreps);
  moinfo.fruocc = init_int_array(nirreps);
  psio_read_entry(CC_INFO, "Frozen Core Orbs Per Irrep",
		  (char *) moinfo.frdocc, sizeof(int)*nirreps);
  psio_read_entry(CC_INFO, "Frozen Virt Orbs Per Irrep",
		  (char *) moinfo.fruocc, sizeof(int)*nirreps);
  
  psio_read_entry(CC_INFO, "No. of Active Orbitals", (char *) &(nactive),
		  sizeof(int)); 

  if(params.ref == 2) {

    moinfo.aoccpi = init_int_array(nirreps);
    moinfo.boccpi = init_int_array(nirreps);
    moinfo.avirtpi = init_int_array(nirreps);
    moinfo.bvirtpi = init_int_array(nirreps);

    psio_read_entry(CC_INFO, "Active Alpha Occ Orbs Per Irrep",
		    (char *) moinfo.aoccpi, sizeof(int)*moinfo.nirreps);
    psio_read_entry(CC_INFO, "Active Beta Occ Orbs Per Irrep",
		    (char *) moinfo.boccpi, sizeof(int)*moinfo.nirreps);
    psio_read_entry(CC_INFO, "Active Alpha Virt Orbs Per Irrep",
		    (char *) moinfo.avirtpi, sizeof(int)*moinfo.nirreps);
    psio_read_entry(CC_INFO, "Active Beta Virt Orbs Per Irrep",
		    (char *) moinfo.bvirtpi, sizeof(int)*moinfo.nirreps);

    moinfo.aocc_sym = init_int_array(nactive);
    moinfo.bocc_sym = init_int_array(nactive);
    moinfo.avir_sym = init_int_array(nactive);
    moinfo.bvir_sym = init_int_array(nactive);

    psio_read_entry(CC_INFO, "Active Alpha Occ Orb Symmetry",
		    (char *) moinfo.aocc_sym, sizeof(int)*nactive);
    psio_read_entry(CC_INFO, "Active Beta Occ Orb Symmetry",
		    (char *) moinfo.bocc_sym, sizeof(int)*nactive);
    psio_read_entry(CC_INFO, "Active Alpha Virt Orb Symmetry",
		    (char *) moinfo.avir_sym, sizeof(int)*nactive);
    psio_read_entry(CC_INFO, "Active Beta Virt Orb Symmetry",
		    (char *) moinfo.bvir_sym, sizeof(int)*nactive);

  }
  else {

    moinfo.occpi = init_int_array(nirreps);
    moinfo.virtpi = init_int_array(nirreps);
    psio_read_entry(CC_INFO, "Active Occ Orbs Per Irrep",
		    (char *) moinfo.occpi, sizeof(int)*nirreps);
    psio_read_entry(CC_INFO, "Active Virt Orbs Per Irrep",
		    (char *) moinfo.virtpi, sizeof(int)*nirreps);

    moinfo.occ_sym = init_int_array(nactive);
    moinfo.vir_sym = init_int_array(nactive);
    psio_read_entry(CC_INFO, "Active Occ Orb Symmetry",
		    (char *) moinfo.occ_sym, sizeof(int)*nactive);
    psio_read_entry(CC_INFO, "Active Virt Orb Symmetry",
		    (char *) moinfo.vir_sym, sizeof(int)*nactive);
  }

  /* Build orbsym array (for AO-basis BT2) */
  moinfo.orbsym = init_int_array(moinfo.nso);
  for(h=0,q=0; h < nirreps; h++)
    for(p=0; p < moinfo.orbspi[h]; p++)
      moinfo.orbsym[q++] = h;

  C = (double ***) malloc(nirreps * sizeof(double **));
  next = PSIO_ZERO;
  for(h=0; h < nirreps; h++) {
    if(moinfo.orbspi[h] && moinfo.virtpi[h]) {
      C[h] = block_matrix(moinfo.orbspi[h],moinfo.virtpi[h]);
      psio_read(CC_INFO, "RHF/ROHF Active Virtual Orbitals", (char *) C[h][0],
		moinfo.orbspi[h]*moinfo.virtpi[h]*sizeof(double), next, &next);
    }
  }
  moinfo.C = C;

  /* Adjust clsdpi array for frozen orbitals */
  for(i=0; i < nirreps; i++)
    moinfo.clsdpi[i] -= moinfo.frdocc[i];

  moinfo.uoccpi = init_int_array(moinfo.nirreps);
  for(i=0; i < nirreps; i++)
    moinfo.uoccpi[i] = moinfo.orbspi[i] - moinfo.clsdpi[i] -
      moinfo.openpi[i] - moinfo.fruocc[i] -
      moinfo.frdocc[i];

  psio_read_entry(CC_INFO, "Reference Energy", (char *) &(moinfo.eref), 
		  sizeof(double));

  fprintf(outfile,"\n\tNuclear Rep. energy (file30)  = %20.15f\n",moinfo.enuc);
  fprintf(outfile,  "\tSCF energy          (file30)  = %20.15f\n",moinfo.escf);
  fprintf(outfile,  "\tReference energy    (file100) = %20.15f\n",moinfo.eref);
}

/* Frees memory allocated in get_moinfo() and dumps out the energy. */
void cleanup(void)
{
  int i, h;

  psio_write_entry(CC_INFO, "CCSD Energy", (char *) &(moinfo.ecc),
		   sizeof(double));

  for(h=0; h < moinfo.nirreps; h++)
    if(moinfo.orbspi[h] && moinfo.virtpi[h]) free_block(moinfo.C[h]);
  free(moinfo.C);

  free(moinfo.orbspi);
  free(moinfo.orbsym);
  free(moinfo.clsdpi);
  free(moinfo.openpi);
  free(moinfo.uoccpi);
  free(moinfo.fruocc);
  free(moinfo.frdocc);
  for(i=0; i < moinfo.nirreps; i++)
    free(moinfo.labels[i]);
  free(moinfo.labels);
  if(params.ref == 2) {
    free(moinfo.aoccpi);
    free(moinfo.boccpi);
    free(moinfo.avirtpi);
    free(moinfo.bvirtpi);
    free(moinfo.aocc_sym);
    free(moinfo.bocc_sym);
    free(moinfo.avir_sym);
    free(moinfo.bvir_sym);
  }
  else {
    free(moinfo.occpi);
    free(moinfo.virtpi);
    free(moinfo.occ_sym);
    free(moinfo.vir_sym);
  }

}

