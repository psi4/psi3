#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<libciomr/libciomr.h>
#include<libpsio/psio.h>
#include<libdpd/dpd.h>
#include<libchkpt/chkpt.h>
#include<libint/libint.h>

#include"moinfo.h"
#include"moinfo_corr.h"

#include"defines.h"
#define EXTERN
#include"global.h"
#include"small_fns.h"

void init_ccinfo()
{
  int *cachefiles, **cachelist;

  cachefiles = init_int_array(PSIO_MAXUNIT);
  /* assuming no caching for DPD here */
  cachelist = init_int_matrix(12,12);

  init_moinfo();
  init_moinfo_corr();

  /*--- Initialize DPD library ---*/
  dpd_init(1, Symmetry.nirreps, UserOptions.memory, 0, cachefiles, cachelist, NULL, 
	   2, MOInfo.clsdpi, MOInfo.mo2symblk_occ, MOInfo.orbspi, MOInfo.mo2symblk);

  free_int_matrix(cachelist,12);
  free(cachefiles);

  return;
}


void cleanup_ccinfo()
{
  dpd_close(1);

  return;
}


