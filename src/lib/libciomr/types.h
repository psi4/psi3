
/* $Id$ */
/* $Log$
 * Revision 1.2  2001/08/08 14:31:25  sherrill
 * removed some slashes in comments to make compiler shut up
 *
 * Revision 1.1.1.1  2000/02/04 22:53:24  evaleev
 * Started PSI 3 repository
 *
 * Revision 2.6  1997/09/14 03:28:55  sherrill
 * Added iosize_() to types.h so the function prototype could be passed to
 * flen.c, which was getting the wrong return type.  Also reformatted some
 * of my old code just a little to make it look a bit nicer.
 *
 * Revision 2.5  1997/09/12  13:53:07  crawdad
 * Changing marco name from ULL to PSI_FPTR.
 *
 * Revision 2.4  1997/08/25  21:50:12  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 2.3  1996/06/18  20:47:48  sherrill
 * Add the whole set of int_array routines to int_array.c (replacing
 * init_int_array.c), add block_matrix.c, and add a new function flen
 * which gets the file size for unit number x.
 *
 * Revision 2.2  1995/04/01  20:53:10  fermann
 * changed bytewise file pointers such as first, last and length to long
 * unsigned ints in order to handle up to 4 gigabyte tmp files (striped into
 * individual pieces of less than 2 gigabytes).  added functions li2sec and
 * sec2li for where they are needed.
 *
 * Revision 2.1  1991/06/15  18:32:30  seidl
 * *** empty log message ***
 * */

#ifndef TYPES_H
#define TYPES_H

#include "iomrparam.h"

typedef
struct {
  int junk;
  } r_async_t;

typedef
struct {
  int junk;
  } s_async_t;

typedef
struct {
  char *path;
#if BUFF
  FILE *stream;
#else
  int stream;
#endif
  } sequential_volume_t;

typedef
struct {
  int n;
  int blocksize;
  int last_ioop;
  int verbose;
  PSI_FPTR next;
  PSI_FPTR previous_size; /* added for error checking with unsigned's */ 
  int unit; /* The unit number. */
  PSI_FPTR incount;  /* The number of bytes read. */
  PSI_FPTR outcount;  /* The number of bytes written. */
  sequential_volume_t v[MAX_VOLUME];
  } sequential_t;

typedef
struct {
  int junk;
  } ram_t;
  

typedef
struct {
  int method;
  union {
    r_async_t *r_async;
    s_async_t *s_async;
    sequential_t *sequential;
    ram_t *ram;
    } ptr;
  } ioFILE_t;


/* From ioopen.c */
void ioinit_();
void ioopen_();
void ioclos_();
void iowrr_();
void iordr_();
void ioabort();
PSI_FPTR iosize_();

/* From sequential.c */
sequential_t *sequential_ioopen();
void sequential_ioclos();
void sequential_iordr();
void sequential_iowrr();
void sequential_iordwrr();
PSI_FPTR sequential_iosize();

/* From r_async.c */
r_async_t *r_async_ioopen();
void r_async_ioclos();
void r_async_iordr();
void r_async_iowrr();
PSI_FPTR r_async_iosize();

/* From s_async.c */
s_async_t *s_async_ioopen();
void s_async_ioclos();
void s_async_iordr();
void s_async_iowrr();
PSI_FPTR s_async_iosize();

/* From ram.c */
ram_t *ram_ioopen();
void ram_ioclos();
void ram_iordr();
void ram_iowrr();
PSI_FPTR ram_iosize();

/* From errors.c */
void no_path_given();
void malloc_check();
void fopen_check();
void fread_error();
void fwrite_error();

#endif
