#include <stdio.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#define EXTERN
#include "globals.h"

void setup(void)
{
  int h, foffset, loffset;

  chkpt_init(PSIO_OPEN_OLD);
  natom = chkpt_rd_natom();
  nmo = chkpt_rd_nmo();
  nso = chkpt_rd_nso();
  nao = chkpt_rd_nao();
  evals = chkpt_rd_evals();
  nirreps = chkpt_rd_nirreps();
  orbspi = chkpt_rd_orbspi();
  clsdpi = chkpt_rd_clsdpi();
  zvals = chkpt_rd_zvals();
  scf = chkpt_rd_scf();
  usotao = chkpt_rd_usotao();
  geom = chkpt_rd_geom();
  chkpt_close();

  ntri = nmo * (nmo + 1)/2;
  ntei = ntri * (ntri + 1)/2;
  noei = nso * (nso + 1)/2;
  noei_ao = nao * (nao + 1)/2;

  uoccpi = init_int_array(nirreps);
  for(h=0; h < nirreps; h++) uoccpi[h] = orbspi[h] - clsdpi[h];

  ndocc = 0; nuocc = 0;
  for(h=0; h < nirreps; h++) { ndocc += clsdpi[h]; nuocc += uoccpi[h]; }

  num_ai = ndocc * nuocc; 
  num_ij = ndocc * (ndocc + 1)/2;

  /* Build the first and last lookup arrays */
  first = init_int_array(nirreps);
  last = init_int_array(nirreps);
  foffset = 0;
  loffset = orbspi[0] - 1;
  first[0] = foffset;
  last[0] = loffset;
  for(h=1; h < nirreps; h++) {
    foffset += orbspi[h-1];
    loffset += orbspi[h];
    first[h] = foffset;
    last[h] = loffset;
  }

  ofirst = init_int_array(nirreps);
  olast = init_int_array(nirreps);
  foffset = 0;
  loffset = clsdpi[0] - 1;
  ofirst[0] = foffset;
  olast[0] = loffset;
  for(h=1; h < nirreps; h++) {
    foffset += orbspi[h-1];
    loffset += uoccpi[h-1] + clsdpi[h];
    ofirst[h] = foffset;
    olast[h] = loffset;
  }

  vfirst = init_int_array(nirreps);
  vlast = init_int_array(nirreps);
  foffset = clsdpi[0];
  loffset = orbspi[0] - 1;
  vfirst[0] = foffset;
  vlast[0] = loffset;
  for(h=1; h < nirreps; h++) {
    foffset += uoccpi[h-1] + clsdpi[h];
    loffset += orbspi[h];
    vfirst[h] = foffset;
    vlast[h] = loffset;
  }
}

void cleanup(void)
{
  free(orbspi);
  free(clsdpi);
  free(uoccpi);
  free(evals);
  free(zvals);
  free(first);
  free(last);
  free(ofirst);
  free(olast);
  free(vfirst);
  free(vlast);
  free(ioff);
  free_block(scf);
  free_block(geom);
}
