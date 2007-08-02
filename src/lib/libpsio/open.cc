/*!
   \file open.cc
   \ingroup (PSIO)
*/

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

using namespace psi;

void
PSIO::open(unsigned int unit, int status)
{
  unsigned int i, j;
  int stream;
  char *name, *path;
  psio_ud *this_unit;

  /* check for too large unit */
  if (unit > PSIO_MAXUNIT)
    psio_error(unit,PSIO_ERROR_MAXUNIT);

  this_unit = &(psio_unit[unit]);

  /* Check to see if this unit is aleady open */
  if(this_unit->vol[0].stream != -1) psio_error(unit,PSIO_ERROR_REOPEN);

  /* Get number of volumes to stripe across */
  this_unit->numvols = get_numvols(unit);
  if(this_unit->numvols > PSIO_MAXVOL) psio_error(unit,PSIO_ERROR_MAXVOL);

  if(!(this_unit->numvols)) this_unit->numvols = 1;

  /* Get the file name prefix */
  get_filename(unit,&name);

  /* Build the name for each volume and open the file */
  for(i=0; i < this_unit->numvols; i++) {
    char* fullpath;
    get_volpath(unit, i, &path);

    if(this_unit->numvols > 1)
      psio_error(unit,PSIO_ERROR_NOVOLPATH);

    fullpath = (char*) malloc( (strlen(path)+strlen(name)+80)*sizeof(char) );
    sprintf(fullpath, "%s%s.%u", path, name, unit);
    this_unit->vol[i].path = strdup(fullpath);
    free(fullpath);

    /* Check if any previously opened volumes have the same path */
    for(j=0; j < i; j++)
      if (!strcmp(this_unit->vol[i].path,this_unit->vol[j].path))
	psio_error(unit,PSIO_ERROR_IDENTVOLPATH);

    /* Now open the volume */
    if(status == PSIO_OPEN_OLD) {
      this_unit->vol[i].stream =
	::open(this_unit->vol[i].path,O_CREAT|O_RDWR,0644);
      if(this_unit->vol[i].stream == -1)
	psio_error(unit,PSIO_ERROR_OPEN);
    }
    else if(status == PSIO_OPEN_NEW) {
      this_unit->vol[i].stream =
	::open(this_unit->vol[i].path,O_CREAT|O_RDWR|O_TRUNC,0644);
      if(this_unit->vol[i].stream == -1)
	psio_error(unit,PSIO_ERROR_OPEN);
    }
    else psio_error(unit,PSIO_ERROR_OSTAT);

    free(path);
  }

  if (status == PSIO_OPEN_OLD) tocread(unit);
  else if (status == PSIO_OPEN_NEW) {
    /* Init the TOC stats and write them to disk */
    this_unit->toclen = 0;
    this_unit->toc = NULL;
    wt_toclen(unit, 0);
  }
  else psio_error(unit,PSIO_ERROR_OSTAT);

  free(name);
}

extern "C" {
  /*!
  ** PSIO_OPEN(): Opens a multivolume PSI direct access file for
  ** reading/writing data.
  **
  **  \param unit   = The PSI unit number used to identify the file to all
  **                  read and write functions.
  **  \param status = Indicates if the file is old (PSIO_OPEN_OLD) or new
  **                  (PSIO_OPEN_NEW). 
  **
  ** \ingroup (PSIO)
  */
  int psio_open(unsigned int unit, int status)
  {
    _default_psio_lib_->open(unit,status);
    return 1;
  }
}

