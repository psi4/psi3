/*!
  \file contr.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <libciomr/libciomr.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_contr(): Reads in the normalized contraction coefficients.
**
**  takes no arguments.
**
**  returns: double *contr Normalized contraction coefficients are
**  returned as an array of doubles. In the checkpoint file they are
**  stored as a matrix MAXANGMOM by the total number of primitives
**  nprim, but each primitive Gaussian contributes to only one shell
**  (and one basis function, of course), so most of these values are
**  zero and not returned.
** \ingroup (CHKPT)
*/

double *chkpt_rd_contr(void)
{
  double *contr;
  double *temp_contr;
  int nprim, i, j, ij = 0;

  nprim = chkpt_rd_nprim();

  temp_contr = init_array(MAXANGMOM*nprim);
  contr = init_array(nprim);

  psio_read_entry(PSIF_CHKPT, "::Contraction coefficients", (char *) temp_contr,
		  MAXANGMOM*nprim*sizeof(double));

/* Picking non-zero coefficients to the "master" array contr */
  for(i=0; i < MAXANGMOM; i++) 
   for(j=0; j < nprim; j++)
    { if (temp_contr[ij] != 0)
         contr[j] = temp_contr[ij];
      ij++;
    }

  free(temp_contr);

  return contr;
}


/*!
** chkpt_wt_contr(): Write out the normalized contraction coefficients.
**
**  \param contr = The array of contraction coefficients.  The
**                 ordering is that given in cints.
**
**  returns: none
**  \ingroup (CHKPT)
*/

void chkpt_wt_contr(double *contr)
{
  int nprim;
  nprim = chkpt_rd_nprim();

  psio_write_entry(PSIF_CHKPT, "::Contraction coefficients", (char *) contr,
		  MAXANGMOM*nprim*sizeof(double));
}
