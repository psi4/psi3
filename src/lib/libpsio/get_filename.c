/*!
   \file get_filename.c
   \ingroup (PSIO)
*/

#include <stdio.h>
#include <stdlib.h>
#include <libipv1/ip_lib.h>
#include "psio.h"

/*!
** PSIO_GET_FILENAME(): Get the filename for filenumber 'unit'
**
** Returns: 
**   0 if a user-specified filename was found
**   1 if the default "psi" is to be used
**
** \ingroup (PSIO)
*/
int psio_get_filename(ULI unit, char *name)
{
  int errcod;
  char ip_token[PSIO_MAXSTR];
  char *gprgid();

  sprintf(ip_token,":%s:FILES:FILE%u:NAME",gprgid(),unit);
  errcod = ip_data(ip_token,"%s",name,0);
  if(errcod == IPE_OK) return(0);

  sprintf(ip_token,":%s:FILES:DEFAULT:NAME",gprgid());
  errcod = ip_data(ip_token,"%s",name,0);
  if(errcod == IPE_OK) return(0);

  sprintf(ip_token,":DEFAULT:FILES:FILE%u:NAME",unit);
  errcod = ip_data(ip_token,"%s",name,0);
  if(errcod == IPE_OK) return(0);

  sprintf(ip_token,":DEFAULT:FILES:DEFAULT:NAME");
  errcod = ip_data(ip_token,"%s",name,0);
  if(errcod == IPE_OK) return(0);

  sprintf(ip_token,":DEFAULT:NAME");
  errcod = ip_data(ip_token,"%s",name,0);
  if(errcod == IPE_OK) return(0);

  /* check the environment */
  if(getenv("PSI_SCRATCH") != NULL) {
    strcpy(name, getenv("PSI_SCRATCH"));
    return(0);
  }

  /* use a default filename */
  strcpy(name, "psi");

  return(1);
}


/*!
** PSIO_GET_FILENAME_DEFAULT(): Get the default filename
*/
int psio_get_filename_default(char *name)
{
  int errcod;
  char ip_token[PSIO_MAXSTR];

  sprintf(ip_token,":DEFAULT:FILES:DEFAULT:NAME");
  errcod = ip_data(ip_token,"%s",name,0);
  if(errcod == IPE_OK) return(0);

  /* check the environment */
  if(getenv("PSI_SCRATCH") != NULL) {
    strcpy(name, getenv("PSI_SCRATCH"));
    return(0);
  }

  /* use a default filename */
  strcpy(name, "psi");

  return(1);
}
