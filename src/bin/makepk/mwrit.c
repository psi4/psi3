
/* $Log$
 * Revision 1.1  2000/02/04 22:51:34  evaleev
 * Initial revision
 *
/* Revision 1.2  1997/08/25 21:52:28  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  22:06:29  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

mwrit(array,address)
   char *array;
   int address;

   {
      int length,iaddress;

      iaddress = block_locs[address];
      length = block_locs[address+nsect/2];

      rwrit(itap40,(char *) array,sizeof(int)*length,iaddress);
      }

mread(array,address)
   char *array;
   int address;

   {
      int length,iaddress;

      iaddress = block_locs[address];
      length = block_locs[address+nsect/2];

      rread(itap40,(char *) array,sizeof(int)*length,iaddress);
      }
