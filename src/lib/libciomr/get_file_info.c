
/* $Log$
 * Revision 1.1  2000/02/04 22:53:18  evaleev
 * Initial revision
 *
/* Revision 2.6  1997/08/25 21:49:48  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 2.5  1997/06/23  12:25:45  crawdad
 * Multiple changes to libciomr: Moved "param.h" to "iomrparam.h" to avoid
 *     conflicts with similarly named system file under linux.  Corrected type
 *    casting in rread(), rwrit(), sread(), and swrit() functions.  Corrected
 *     unclosed tmpdisks.dat file in sequential.c.  Corrected block_matrix() to
 *    avoid malloc'ing zero-length arrays.
 *
 * -Daniel
 *
 * Revision 2.4  1996/07/02  21:06:27  sherrill
 * Completed the changes needed to use unit numbers greater than 99; changed
 * a hardwired 99 in init_ptrs.c to MAX_UNIT and increased a "unit" string
 * from 3 chars to 4.  Also removed a compiler warning in sequential.c by
 * casting ud to (char *) for malloc_check().
 *
 * Revision 2.3  1991/09/18  20:47:14  seidl
 * dec changes
 *
 * Revision 2.2  1991/08/21  05:41:35  psi
 * declare gprgid
 *
 * Revision 2.1  1991/06/15  18:28:54  seidl
 * initial revision
 * */

static char *rcsid = "$Id$";


#include "iomrparam.h"
#include "includes.h"
#include <ip_libv1.h>

int get_file_info(token,format,val)
  char *token,*format;
#ifdef DEC
  char *val;
#else
  void *val;
#endif
{
  int i,errcod;
  char *prog,unit[4],*junk;
  char ip_token[MAX_STRING];
  char *gprgid();

  prog = gprgid();
  junk = strchr(token,':');
  junk++;
  i=0;
  while(*junk != ':') {
    unit[i] = *junk;
    junk++;
    i++;
    }
  junk++;
  unit[i]='\0';

  sprintf(ip_token,":%s:FILES:FILE%s:%s",prog,unit,junk);
  errcod = ip_data(ip_token,format,val,0);
  if(errcod == IPE_OK) return(0);

  sprintf(ip_token,":%s:FILES:DEFAULT:%s",prog,junk);
  errcod = ip_data(ip_token,format,val,0);
  if(errcod == IPE_OK) return(0);

  sprintf(ip_token,":DEFAULT:FILES:FILE%s:%s",unit,junk);
  errcod = ip_data(ip_token,format,val,0);
  if(errcod == IPE_OK) return(0);

  sprintf(ip_token,":DEFAULT:FILES:DEFAULT:%s",junk);
  errcod = ip_data(ip_token,format,val,0);
  if(errcod == IPE_OK) return(0);

  return(-1);
  }
