/*!
   \file get_numvols.c
   \ingroup (PSIO)
*/
 
#include <stdio.h>
#include <string.h>
#include <libipv1/ip_lib.h>
#include "psio.h"

/*!
** PSIO_GET_NUMVOLS(): Get the number of volumes that file number 'unit'
** is split across.
**
** \ingroup (PSIO)
*/
ULI psio_get_numvols(ULI unit)
{
  ULI num;
  int errcod;
  char ip_token[PSIO_MAXSTR];
  char *gprgid();

  sprintf(ip_token,":%s:FILES:FILE%u:NVOLUME",gprgid(),unit);
  errcod = ip_data(ip_token,"%u",&num,0);
  if(errcod == IPE_OK) return(num);

  sprintf(ip_token,":%s:FILES:DEFAULT:NVOLUME",gprgid());
  errcod = ip_data(ip_token,"%u",&num,0);
  if(errcod == IPE_OK) return(num);

  sprintf(ip_token,":PSI:FILES:FILE%u:NVOLUME",unit);
  errcod = ip_data(ip_token,"%u",&num,0);
  if(errcod == IPE_OK) return(num);

  sprintf(ip_token,":PSI:FILES:DEFAULT:NVOLUME");
  errcod = ip_data(ip_token,"%u",&num,0);
  if(errcod == IPE_OK) return(num);

  sprintf(ip_token,":DEFAULT:FILES:FILE%u:NVOLUME",unit);
  errcod = ip_data(ip_token,"%u",&num,0);
  if(errcod == IPE_OK) return(num);

  sprintf(ip_token,":DEFAULT:FILES:DEFAULT:NVOLUME");
  errcod = ip_data(ip_token,"%u",&num,0);
  if(errcod == IPE_OK) return(num);

  /* default to one volume */
  return(1);
}


/*!
** PSIO_GET_NUMVOLS_DEFAULT(): Get the number of volumes that file 
** number 'unit' is split across.
**
** \ingroup (PSIO)
*/
ULI psio_get_numvols_default(void)
{
  ULI num;
  int errcod;
  char ip_token[PSIO_MAXSTR];

  sprintf(ip_token,":PSI:FILES:DEFAULT:NVOLUME");
  errcod = ip_data(ip_token,"%u",&num,0);
  if(errcod == IPE_OK) return(num);

  sprintf(ip_token,":DEFAULT:FILES:DEFAULT:NVOLUME");
  errcod = ip_data(ip_token,"%u",&num,0);
  if(errcod == IPE_OK) return(num);

  /* default to one volume */
  return(1);
}
