/*!
   \file get_volpath.c
   \ingroup (PSIO)
*/

#include <stdio.h>
#include <string.h>
#include <libpsio/psio.h>

static const char* default_path = "/tmp/";
extern char *gprgid();

/*
** PSIO_GET_VOLPATH(): Copy the path to a given volume for file number
** 'unit' into 'path' (storage is malloc'ed).
**
** \ingroup (PSIO)
*/
int psio_get_volpath(unsigned int unit, unsigned int volume, char **path)
{
  const char* kval;
  char volumeX[20];
  sprintf(volumeX,"VOLUME%u",volume+1);

  kval = psio_get_filescfg_kwd(gprgid(),volumeX,unit);
  if(kval != 0) { *path = strdup(kval); return(0); }
  kval = psio_get_filescfg_kwd(gprgid(),volumeX,-1);
  if(kval != 0) { *path = strdup(kval); return(0); }
  kval = psio_get_filescfg_kwd("PSI",volumeX,unit);
  if(kval != 0) { *path = strdup(kval); return(0); }
  kval = psio_get_filescfg_kwd("PSI",volumeX,-1);
  if(kval != 0) { *path = strdup(kval); return(0); }
  kval = psio_get_filescfg_kwd("DEFAULT",volumeX,unit);
  if(kval != 0) { *path = strdup(kval); return(0); }
  kval = psio_get_filescfg_kwd("DEFAULT",volumeX,-1);
  if(kval != 0) { *path = strdup(kval); return(0); }

  *path = strdup(default_path);
  return(1);
}


/*
** PSIO_GET_VOLPATH_DEFAULT(): Get the default path for the nth volume
** of any file.
**
** \ingroup (PSIO)
*/
int psio_get_volpath_default(unsigned int volume, char **path)
{
  const char* kval;
  char volumeX[20];
  sprintf(volumeX,"VOLUME%u",volume+1);

  kval = psio_get_filescfg_kwd("PSI",volumeX,-1);
  if(kval != 0) { *path = strdup(kval); return(0); }
  kval = psio_get_filescfg_kwd("DEFAULT",volumeX,-1);
  if(kval != 0) { *path = strdup(kval); return(0); }

  *path = strdup(default_path);
  return(1);
}
