
/* $Log$
 * Revision 1.1  2000/02/04 22:53:20  evaleev
 * Initial revision
 *
/* Revision 2.7  1999/11/01 20:10:56  evaleev
/* Added explicit extern declarations of functions within the library.
/*
/* Revision 2.6  1997/09/12 13:52:52  crawdad
/* Changing marco name from ULL to PSI_FPTR.
/*
 * Revision 2.5  1997/08/25  21:49:56  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 2.4  1997/06/23  12:25:51  crawdad
 * Multiple changes to libciomr: Moved "param.h" to "iomrparam.h" to avoid
 *     conflicts with similarly named system file under linux.  Corrected type
 *    casting in rread(), rwrit(), sread(), and swrit() functions.  Corrected
 *     unclosed tmpdisks.dat file in sequential.c.  Corrected block_matrix() to
 *    avoid malloc'ing zero-length arrays.
 *
 * -Daniel
 *
 * Revision 2.3  1996/06/18  20:47:40  sherrill
 * Add the whole set of int_array routines to int_array.c (replacing
 * init_int_array.c), add block_matrix.c, and add a new function flen
 * which gets the file size for unit number x.
 *
 * Revision 2.2  1995/04/01  20:52:55  fermann
 * changed bytewise file pointers such as first, last and length to long
 * unsigned ints in order to handle up to 4 gigabyte tmp files (striped into
 * individual pieces of less than 2 gigabytes).  added functions li2sec and
 * sec2li for where they are needed.
 *
 * Revision 2.1  1991/06/15  18:29:25  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id$";


#define MAIN

#include "iomrparam.h"
#include "includes.h"
#include "types.h"

#ifdef DEC
extern int get_param(char *, char *, char *);
extern int get_file_info(char *, char *, char *);
#else
extern int get_param(char *, char *, void *);
extern int get_file_info(char *, char *, void *);
#endif
extern sequential_t *sequential_ioopen(char *, int);
extern void sequential_ioclos(sequential_t *, int);
extern void sequential_iordr(sequential_t *, char *, PSI_FPTR, int);
extern void sequential_iowrr(sequential_t *, char *, PSI_FPTR, int);
extern void sequential_iordwrr(char *, int, sequential_t *, char *, PSI_FPTR, int);
extern PSI_FPTR sequential_iosize(sequential_t *);
extern r_async_t *r_async_ioopen(char *, int);
extern void r_async_ioclos(r_async_t *, int);
extern void r_async_iordr(r_async_t *, char *, PSI_FPTR, int);
extern void r_async_iowrr(r_async_t *, char *, PSI_FPTR, int);
extern PSI_FPTR r_async_iosize(r_async_t *);
extern s_async_t *s_async_ioopen(char *, int);
extern void s_async_ioclos(s_async_t *, int);
extern void s_async_iordr(s_async_t *, char *, PSI_FPTR, int);
extern void s_async_iowrr(s_async_t *, char *, PSI_FPTR, int);
extern PSI_FPTR s_async_iosize(s_async_t *);
extern ram_t *ram_ioopen(char *, int);
extern void ram_ioclos(ram_t *, int);
extern void ram_iordr(ram_t *, char *, PSI_FPTR, int);
extern void ram_iowrr(ram_t *, char *, PSI_FPTR, int);
extern PSI_FPTR ram_iosize(ram_t *);

ioFILE_t ud[MAX_UNIT];
unsigned int current_unit;

void
ioinit_()
{
  int i;

  for (i=0; i<MAX_UNIT; i++) {
    ud[i].method = PSIIO_UNOPENED;
    ud[i].ptr.sequential = NULL;
    }
  }

void
ioopen_(unit)
int *unit;
{
  char param[MAX_STRING];
  char methodparam[MAX_STRING];
  char method[MAX_STRING]; 
  char unitch[MAX_STRING];

  current_unit = *unit;

  sprintf(unitch,"%d",*unit);

  strcpy(param,"FILES:");
  strcat(param,unitch);
  strcat(param,":");

  strcpy(methodparam,param);

  if(oldstyleinput()) {
    strcat(methodparam,"method");
    if (get_param(methodparam,"%s",method) != 0) {
      strcpy(method,"sequential");
      }
    }
  else {
    strcat(methodparam,"METHOD");
    if (get_file_info(methodparam,"%s",method) != 0) {
      strcpy(method,"sequential");
      }
    }

  if (!strcmp(method,"sequential")) {
    ud[*unit].method = PSIIO_SEQUENTIAL;
    ud[*unit].ptr.sequential = sequential_ioopen(param,*unit);
    }
  else if (!strcmp(method,"r_async")) {
    ud[*unit].method = PSIIO_R_ASYNC;
    ud[*unit].ptr.r_async = r_async_ioopen(param,*unit);
    }
  else if (!strcmp(method,"s_async")) {
    ud[*unit].method = PSIIO_S_ASYNC;
    ud[*unit].ptr.s_async = s_async_ioopen(param,*unit);
    }
  else if (!strcmp(method,"ram")) {
    ud[*unit].method = PSIIO_RAM;
    ud[*unit].ptr.ram = ram_ioopen(param,*unit);
    }
  else {
    fprintf(stderr,"ioopen_: invalid method (%s) for unit %d\n",method,*unit);
    ioabort();
    }
  }

void
ioclos_(unit,status)
int *unit, *status;
{

  current_unit = *unit;

  if (ud[*unit].method == PSIIO_SEQUENTIAL) {
    ud[*unit].method = PSIIO_UNOPENED;
    sequential_ioclos(ud[*unit].ptr.sequential, *status);
    }
  else if (ud[*unit].method == PSIIO_R_ASYNC) {
    ud[*unit].method = PSIIO_UNOPENED;
    r_async_ioclos(ud[*unit].ptr.r_async, *status);
    }
  else if (ud[*unit].method == PSIIO_S_ASYNC) {
    ud[*unit].method = PSIIO_UNOPENED;
    s_async_ioclos(ud[*unit].ptr.s_async, *status);
    }
  else if (ud[*unit].method == PSIIO_RAM) {
    ud[*unit].method = PSIIO_UNOPENED;
    ram_ioclos(ud[*unit].ptr.ram, *status);
    }
  else {
    fprintf(stderr,"ioclos_: invalid method (%d) for unit %d\n",
             ud[*unit].method,*unit);
    ioabort();
    }
  free(ud[*unit].ptr.sequential);
  ud[*unit].ptr.sequential = NULL;
  ud[*unit].method = PSIIO_UNOPENED;
  }

void
iowrr_(unit,buffer,first,length)
int *unit,*length;
PSI_FPTR *first;
char *buffer;
{

  current_unit = *unit;

  if (*first + *length < *first) {
    fprintf(stderr,"iowrr_: *unit=%d, *first=%lu, *length=%lu\n",
            *unit, *first, *length);
    ioabort();
    }

  if (ud[*unit].method == PSIIO_SEQUENTIAL) {
    sequential_iowrr(ud[*unit].ptr.sequential,buffer,*first,*length);
    }
  else if (ud[*unit].method == PSIIO_R_ASYNC) {
    r_async_iowrr(ud[*unit].ptr.r_async,buffer,*first,*length);
    }
  else if (ud[*unit].method == PSIIO_S_ASYNC) {
    s_async_iowrr(ud[*unit].ptr.s_async,buffer,*first,*length);
    }
  else if (ud[*unit].method == PSIIO_RAM) {
    ram_iowrr(ud[*unit].ptr.ram,buffer,*first,*length);
    }
  else {
    fprintf(stderr,"iowrr_: invalid method (%d) for unit %d\n",
            ud[*unit].method,*unit);
    ioabort();
    }
  }

void
iordr_(unit,buffer,first,length)
int *unit, *length;
PSI_FPTR *first;
char *buffer;
{

  current_unit = *unit;

  if (*first + *length < *first) {
    fprintf(stderr,"iordr_: *unit=%d, *first=%lu, *length=%lu\n",
            *unit, *first, *length);
    ioabort();
    }

  if (ud[*unit].method == PSIIO_SEQUENTIAL) {
    sequential_iordr(ud[*unit].ptr.sequential,buffer,*first,*length);
    }
  else if (ud[*unit].method == PSIIO_R_ASYNC) {
    r_async_iordr(ud[*unit].ptr.r_async,buffer,*first,*length);
    }
  else if (ud[*unit].method == PSIIO_S_ASYNC) {
    s_async_iordr(ud[*unit].ptr.s_async,buffer,*first,*length);
    }
  else if (ud[*unit].method == PSIIO_RAM) {
    ram_iordr(ud[*unit].ptr.ram,buffer,*first,*length);
    }
  else {
    fprintf(stderr,"iordr_: invalid method (%d) for unit %d\n",
            ud[*unit].method,*unit);
    ioabort();
    }

  }

void
ioabort()
{
  fprintf(stderr,"ioabort: current unit = %d\n",current_unit);
  exit(1);
  }

PSI_FPTR
iosize_(unit)
int *unit;
{

   PSI_FPTR fsize=0;

   current_unit = *unit;

   if (ud[*unit].method == PSIIO_UNOPENED) {
      fprintf(stderr,"iosize_: Unit %d not opened\n", current_unit);
      return(0);
      }

   else if (ud[*unit].method == PSIIO_SEQUENTIAL) {
      fsize = sequential_iosize(ud[*unit].ptr.sequential); 
      }

   else if (ud[*unit].method == PSIIO_R_ASYNC) {
      fsize = r_async_iosize(ud[*unit].ptr.r_async); 
      }

   else if (ud[*unit].method == PSIIO_S_ASYNC) {
      fsize = s_async_iosize(ud[*unit].ptr.s_async); 
      }

   else if (ud[*unit].method == PSIIO_RAM) {
      fsize = ram_iosize(ud[*unit].ptr.ram); 
      }

   else {
      fprintf(stderr,"iosize_: Invalid method for unit %d\n", current_unit);
      return(0);
      }

   return(fsize);
}
      
