/*!
  \file zmat.c
  \ingroup (CHKPT)
*/

#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_zmat():  Reads in the z_matrix.
**
**   takes no arguments.
**
**   returns: z_geom = An array natom long which contains 
**     a z_entry struct for each atom
** 
** \ingroup (CHKPT)
*/

struct z_entry *chkpt_rd_zmat(void)
{
  int nallatom;
  struct z_entry *z_geom;
  char *keyword;
  keyword = chkpt_build_keyword("Z-matrix");

  nallatom = chkpt_rd_nallatom();
  z_geom = (struct z_entry *) malloc(nallatom*(sizeof(struct z_entry)));

  psio_read_entry(PSIF_CHKPT, keyword, (char *) z_geom, 
    sizeof(struct z_entry)*nallatom);

  free(keyword);
  return z_geom;
}


/*!
** chkpt_wt_zmat():  Writes out the z_matrix.
**
**  \param z_geom = An array natom long which contains 
**     a z_entry struct for each atom
**
** returns: none
**
** \ingroup (CHKPT)
*/

void chkpt_wt_zmat(struct z_entry *z_geom)
{
  int nallatom;
  char *keyword;
  keyword = chkpt_build_keyword("Z-matrix");

  nallatom = chkpt_rd_nallatom();

  psio_write_entry(PSIF_CHKPT, keyword, (char *) z_geom, 
    sizeof(struct z_entry)*nallatom);

  free(keyword);
}
