/*!
  \file psi_stop
  \ingroup (CIOMR)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <psifiles.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>

extern FILE *infile, *outfile;
extern char *psi_file_prefix;

/*!
** psi_stop()
**
** This function closes input and output files and deinitializes Input Parsing library.
**
** Arguments: none
**
** Returns: one of standard PSI error codes
** \ingroup (CIOMR)
*/

int psi_stop(void)
{
  ip_done();
  free(psi_file_prefix);
  fclose(outfile);
  fclose(infile);

  return(PSI_RETURN_SUCCESS);
}
