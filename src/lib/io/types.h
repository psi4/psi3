#include "param.h"

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
  FILE *stream;
  } sequential_volume_t;

typedef
struct {
  int n;
  int blocksize;
  int last_ioop;
  int verbose;
  PSI_FPTR next;
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
void io_q_abort();

/* From sequential.c */
sequential_t *sequential_ioopen();
void sequential_ioclos();
void sequential_iordr();
void sequential_iowrr();
void sequential_iordwrr();

/* From r_async.c */
r_async_t *r_async_ioopen();
void r_async_ioclos();
void r_async_iordr();
void r_async_iowrr();

/* From s_async.c */
s_async_t *s_async_ioopen();
void s_async_ioclos();
void s_async_iordr();
void s_async_iowrr();

/* From ram.c */
ram_t *ram_ioopen();
void ram_ioclos();
void ram_iordr();
void ram_iowrr();

/* From errors.c */
void no_path_given();
void malloc_check();
void fopen_check();
void fread_error();
void fwrite_error();
