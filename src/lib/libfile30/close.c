#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_close()  closes up file30, frees file30 struct.
** 
**  arguments: none, but file30_init must already have been called for 
**    this to work.  
**
**  returns: zero.  Perhaps this, too, will change one day.
*/


int file30_close(void)
{
  /* Free the arrays in info30_ */
  free(info30_.mpoint);
  free(info30_.mconst);
  free(info30_.mcalcs);
  info30_.mpoint = NULL;
  info30_.mconst = NULL;
  info30_.mcalcs = NULL;

  /* Close up the file */
  rclose(info30_.filenum,3);

  return 0;
}
