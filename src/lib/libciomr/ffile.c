/* $Log$
 * Revision 1.3  2002/06/01 18:23:54  sherrill
 * Upgrade doxygen documentation
 *
/* Revision 1.2  2002/04/18 21:47:35  sherrill
/* Here's some changes to document via doxygen and upgrade to ANSI C
/* instead of K&R declarations.
/*
/* Revision 1.1.1.1  2000/02/04 22:53:18  evaleev
/* Started PSI 3 repository
/*
/* Revision 2.3  1997/08/25 21:49:46  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 2.2  1997/06/23  12:25:44  crawdad
 * Multiple changes to libciomr: Moved "param.h" to "iomrparam.h" to avoid
 *     conflicts with similarly named system file under linux.  Corrected type
 *    casting in rread(), rwrit(), sread(), and swrit() functions.  Corrected
 *     unclosed tmpdisks.dat file in sequential.c.  Corrected block_matrix() to
 *    avoid malloc'ing zero-length arrays.
 *
 * -Daniel
 *
 * Revision 2.1  1991/06/15  18:28:52  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id$";

/*!
** \file ffile.c
** \ingroup (CIOMR)
*/ 

#include "iomrparam.h"
#include "includes.h"

/*!
** ffile(): Open an ASCII file for reading/writing.  Returns a pointer
** to the new file in argument unit.
**
** \param unit = pointer to hold a FILE pointer
** \param name = filename to open
** \param code = 0 (write), 1 (write/append), 2 (read)
** \ingroup (CIOMR)
*/
void ffile(FILE **unit,char *name,int code)
  {
      if (name[0] == ' ') {
         fprintf(stderr,"ffile no longer has a default filename\n");
         exit(2);
         }

      switch (code) {
         case 0:
            if ((*unit = fopen(name,"w+"))==NULL){
               fprintf(stderr,"error in ffile: cannot open file %s\n",name);
               exit(2);
               }
            break;
         case 1:
            if ((*unit = fopen(name,"a+"))==NULL){
               fprintf(stderr,"error in ffile: cannot open file %s\n",name);
               exit(2);
               }
            break;
         case 2:
            if ((*unit = fopen(name,"r+"))==NULL){
               fprintf(stderr,"error in ffile: cannot open file %s\n",name);
               exit(2);
               }
            break;
         default:
            fprintf(stderr,"error in ffile: invalid code %d\n",code);
         }
      }

