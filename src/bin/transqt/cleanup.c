#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>
#include "MOInfo.h"
#include "Params.h"
#include "globals.h"

void destruct_evects(int nirreps, double ***evects);


void cleanup(void)
{
  int i;

  /* Free moinfo Arrays */
  free(moinfo.sopi);
  free(moinfo.orbspi);
  free(moinfo.clsdpi);
  free(moinfo.openpi);
  free(moinfo.virtpi);
  free(moinfo.sosym);
  free(moinfo.orbsym);
  free(moinfo.order);
  free(moinfo.corr2pitz);
  free(moinfo.fruocc);
  if (params.backtr) free(moinfo.corr2pitz_nofzv);
  free(moinfo.frdocc);
  free(moinfo.first_so);
  free(moinfo.last_so);
  free(moinfo.first);
  free(moinfo.last);
  free(moinfo.fstact);
  free(moinfo.lstact);
  for(i=0; i < moinfo.nirreps; i++)
      free(moinfo.labels[i]);
  free(moinfo.labels);
  destruct_evects(params.backtr ? moinfo.backtr_nirreps : moinfo.nirreps, 
                  moinfo.evects);
  free(moinfo.active);
  free_block(moinfo.scf_vector);
  free(moinfo.evals);
  free(moinfo.oe_ints);
  free(moinfo.fzc_operator);
  free(moinfo.S);
  if (params.reorder) free(params.moorder);
  free(moinfo.backtr_mo_first);
  free(moinfo.backtr_mo_last);
  free(moinfo.backtr_mo_fstact);
  free(moinfo.backtr_mo_lstact);
  free(moinfo.backtr_mo_orbspi);
  free(moinfo.backtr_mo_active);
  free(moinfo.backtr_ao_first);
  free(moinfo.backtr_ao_last);
  free(moinfo.backtr_ao_orbspi);
  free(moinfo.backtr_ao_orbsym);

  /* Free ioff Array */
  free(ioff);

  /* Free params Arrays */
  free(params.wfn);
}
