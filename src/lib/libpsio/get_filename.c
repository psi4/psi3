#include <stdio.h>
#include <libipv1/ip_lib.h>
#include "psio.h"

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

  return(1);
}


int psio_get_filename_default(char *name)
{
  int errcod;
  char ip_token[PSIO_MAXSTR];

  sprintf(ip_token,":DEFAULT:FILES:DEFAULT:NAME");
  errcod = ip_data(ip_token,"%s",name,0);
  if(errcod == IPE_OK) return(0);

  return(1);
}
