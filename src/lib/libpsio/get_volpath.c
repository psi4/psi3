#include <stdio.h>
#include <libipv1/ip_lib.h>
#include "psio.h"

int psio_get_volpath(ULI unit, ULI volume, char *path)
{
  int errcod;
  char ip_token[PSIO_MAXSTR];
  char *gprgid();

  sprintf(ip_token,":%s:FILES:FILE%u:VOLUME%u",gprgid(),unit,volume+1);
  errcod = ip_data(ip_token,"%s",path,0);
  if(errcod == IPE_OK) return(0);

  sprintf(ip_token,":%s:FILES:DEFAULT:VOLUME%u",gprgid(),volume+1);
  errcod = ip_data(ip_token,"%s",path,0);
  if(errcod == IPE_OK) return(0);

  sprintf(ip_token,":DEFAULT:FILES:FILE%u:VOLUME%u",unit,volume+1);
  errcod = ip_data(ip_token,"%s",path,0);
  if(errcod == IPE_OK) return(0);

  sprintf(ip_token,":DEFAULT:FILES:DEFAULT:VOLUME%u",volume+1);
  errcod = ip_data(ip_token,"%s",path,0);
  if(errcod == IPE_OK) return(0);

  return(1);
}


int psio_get_volpath_default(ULI volume, char *path)
{
  int errcod;
  char ip_token[PSIO_MAXSTR];

  sprintf(ip_token,":DEFAULT:FILES:DEFAULT:VOLUME%u",volume+1);
  errcod = ip_data(ip_token,"%s",path,0);
  if(errcod == IPE_OK) return(0);

  return(1);
}
