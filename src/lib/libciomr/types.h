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
