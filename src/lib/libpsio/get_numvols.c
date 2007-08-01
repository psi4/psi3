/*!
   \file get_numvols.c
   \ingroup (PSIO)
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <libpsio/psio.h>

extern char *gprgid();

/*!
** PSIO_GET_NUMVOLS(): Get the number of volumes that file number 'unit'
** is split across.
**
** \ingroup (PSIO)
*/
unsigned int psio_get_numvols(unsigned int unit)
{
  const char* charnum;

  charnum = (char*) psio_get_filescfg_kwd(gprgid(),"NVOLUME",unit);
  if(charnum != 0) return((unsigned int)atoi(charnum));
  charnum = (char*) psio_get_filescfg_kwd(gprgid(),"NVOLUME",-1);
  if(charnum != 0) return((unsigned int)atoi(charnum));
  charnum = (char*) psio_get_filescfg_kwd("PSI","NVOLUME",unit);
  if(charnum != 0) return((unsigned int)atoi(charnum));
  charnum = (char*) psio_get_filescfg_kwd("PSI","NVOLUME",-1);
  if(charnum != 0) return((unsigned int)atoi(charnum));
  charnum = (char*) psio_get_filescfg_kwd("DEFAULT","NVOLUME",unit);
  if(charnum != 0) return((unsigned int)atoi(charnum));
  charnum = (char*) psio_get_filescfg_kwd("DEFAULT","NVOLUME",-1);
  if(charnum != 0) return((unsigned int)atoi(charnum));

  /* default to one volume */
  return(1);
}


/*!
** PSIO_GET_NUMVOLS_DEFAULT(): Get the number of volumes that file 
** number 'unit' is split across.
**
** \ingroup (PSIO)
*/
unsigned int psio_get_numvols_default(void)
{
  const char* charnum;

  charnum = (char*) psio_get_filescfg_kwd("PSI","NVOLUME",-1);
  if(charnum != 0) return((unsigned int)atoi(charnum));
  charnum = (char*) psio_get_filescfg_kwd("DEFAULT","NVOLUME",-1);
  if(charnum != 0) return((unsigned int)atoi(charnum));

  /* default to one volume */
  return(1);
}
