/*!
  \file rd_shell_per_am.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_rd_shells_per_am(): Reads in the numbers of shells of each angular
**                            momentum.
**
**  takes no arguments.
**
**  returns: int *shells_per_am
*/


int *file30_rd_shells_per_am(void)
{
  int *spam;
  int max_am;
  PSI_FPTR spam_ptr, junk;

  max_am = file30_rd_max_am();
  spam_ptr = (PSI_FPTR) (info30_.mpoint[47] - 1)*sizeof(int);
  spam = init_int_array(max_am+1);

  wreadw(info30_.filenum, (char *) spam, (int) (max_am+1)*sizeof(int),
	 spam_ptr, &junk);

  return spam;
}
