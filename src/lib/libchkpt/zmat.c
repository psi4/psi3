/*!
  \file zmat.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <libfile30/file30.h>
#include <psifiles.h>
#include <stdlib.h>
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
  int nentry;
  struct z_entry *z_geom;

  nentry = chkpt_rd_nentry();

  z_geom = (struct z_entry *) malloc(nentry*(sizeof(struct z_entry)));

  psio_read_entry(PSIF_CHKPT, "::Z-matrix", (char *) z_geom, 
                  sizeof(struct z_entry)*nentry);

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
  int nentry;

  nentry = chkpt_rd_nentry();

  psio_write_entry(PSIF_CHKPT, "::Z-matrix", (char *) z_geom, 
                   sizeof(struct z_entry)*nentry);
}
