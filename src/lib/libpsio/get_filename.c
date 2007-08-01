/*!
   \file get_filename.c
   \ingroup (PSIO)
*/

#include <stdio.h>
#include <string.h>
#include <libpsio/psio.h>

extern char *psi_file_prefix;
extern char *gprgid();

/*!
** PSIO_GET_FILENAME(): Copy the filename for filenumber 'unit' into 'name' (it is malloc'ed)
**
** Returns: 
**   0 if a user-specified filename was found
**   1 if the global default will be used
**
** \ingroup (PSIO)
*/
int psio_get_filename(unsigned int unit, char **name)
{
  const char* kval;
  kval = psio_get_filescfg_kwd(gprgid(),"NAME",unit);
  if(kval != 0) { *name = strdup(kval); return(0); }
  kval = psio_get_filescfg_kwd(gprgid(),"NAME",-1);
  if(kval != 0) { *name = strdup(kval); return(0); }
  kval = psio_get_filescfg_kwd("PSI","NAME",unit);
  if(kval != 0) { *name = strdup(kval); return(0); }
  kval = psio_get_filescfg_kwd("PSI","NAME",-1);
  if(kval != 0) { *name = strdup(kval); return(0); }
  kval = psio_get_filescfg_kwd("DEFAULT","NAME",unit);
  if(kval != 0) { *name = strdup(kval); return(0); }
  kval = psio_get_filescfg_kwd("DEFAULT","NAME",-1);
  if(kval != 0) { *name = strdup(kval); return(0); }

  *name = strdup(psi_file_prefix);
  return(1);
}


/*!
** PSIO_GET_FILENAME_DEFAULT(): Get the default filename
*/
int psio_get_filename_default(char **name)
{
  *name = strdup(psi_file_prefix);
  return(1);
}
