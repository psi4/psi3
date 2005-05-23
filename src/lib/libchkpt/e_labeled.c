/*!
  \file e_labeled.c
  \ingroup (CHKPT)
*/

#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_e_labeled(): Reads in an energy with a given label
**
**  arguments: 
**   \param char * label
**
**  returns: double E, the energy
**  \ingroup (CHKPT)
*/

double chkpt_rd_e_labeled(char *label)
{
  char *s;
  double E;

  s = (char *) malloc (strlen(label)+3);
  strcpy(s,"::");
  strcat(s,label);
  /* printf("chkpt_rd_e_labeled using label %s\n",s); */

  psio_read_entry(PSIF_CHKPT, s, (char *) &E, sizeof(double));

  free(s);
  return E;
}

/*!
** chkpt_wt_e_labeled(): Write an energy along with a label
**
**  arguments: 
**   \param char *label, the label
**   \param double E, the energy
**
**  returns: none
**  \ingroup (CHKPT)
*/

void chkpt_wt_e_labeled(char *label, double E)
{
  char *s;

  s = (char *) malloc (strlen(label)+3);
  strcpy(s,"::");
  strcat(s,label);
  /* printf("chkpt_wt_e_labeled using label %s\n",s); */

  psio_write_entry(PSIF_CHKPT, s, (char *) &E, sizeof(double));
  free(s);
}

