
#define MAIN

#include "includes.h"
#include "param.h"
#include "types.h"

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
C_IOOPEN (unit)
int *unit;
{
  char param[MAX_STRING];
  char methodparam[MAX_STRING];
  char method[MAX_STRING]; 
  char unitch[MAX_STRING];

  current_unit = *unit;

  sprintf(unitch,"%d",*unit);

  if (oldstyleinput()) {
    strcpy(param,"FILES:");
    strcat(param,unitch);
    strcat(param,":");

    strcpy(methodparam,param);
    strcat(methodparam,"method");

    if (get_param(methodparam,"%s",method) != 0) {
      strcpy(method,"sequential");
      }
    }
  else {
    if (!read_files_string(*unit,"METHOD",method)) {
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
  else if (!strcmp(method,"fortran")) {
    ud[*unit].method = PSIIO_FORTRAN;
    C_IOFO (unit);
    }
  else {
    fprintf(stderr,"ioopen_: invalid method (%s) for unit %d\n",method,*unit);
    io_q_abort();
    }
  }

void
C_IOCLOS (unit,status)
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
  else if (ud[*unit].method == PSIIO_FORTRAN) {
    ud[*unit].method = PSIIO_UNOPENED;
    C_IOFC (unit, status);
    }
  else {
    fprintf(stderr,"ioclos_: invalid method (%d) for unit %d\n",
             ud[*unit].method,*unit);
    io_q_abort();
    }
  free(ud[*unit].ptr.sequential);
  ud[*unit].ptr.sequential = NULL;
  ud[*unit].method = PSIIO_UNOPENED;
  }

void
C_IOWRR (unit,buffer,ifirst,ilength)
int *unit;
unsigned long int *ifirst,*ilength;
char *buffer;
{
  unsigned long int length;
  PSI_FPTR first;

  first = (PSI_FPTR) (*ifirst-1)*sizeof(int);
  length = (*ilength)*sizeof(int);

  current_unit = *unit;

  if (first + length < first) {
    fprintf(stderr,"iowrr_: *unit=%d, first=%lu, length=%lu\n",
            *unit, first, length);
    io_q_abort();
    }

  if (ud[*unit].method == PSIIO_SEQUENTIAL) {
    sequential_iowrr(ud[*unit].ptr.sequential,buffer,first,length);
    }
  else if (ud[*unit].method == PSIIO_R_ASYNC) {
    r_async_iowrr(ud[*unit].ptr.r_async,buffer,first,length);
    }
  else if (ud[*unit].method == PSIIO_S_ASYNC) {
    s_async_iowrr(ud[*unit].ptr.s_async,buffer,first,length);
    }
  else if (ud[*unit].method == PSIIO_RAM) {
    ram_iowrr(ud[*unit].ptr.ram,buffer,first,length);
    }
  else if (ud[*unit].method == PSIIO_FORTRAN) {
    C_IOFW (unit,buffer,&first,&length);
    }
  else {
    fprintf(stderr,"iowrr_: invalid method (%d) for unit %d\n",
            ud[*unit].method,*unit);
    io_q_abort();
    }
  }

void
C_IORDR (unit,buffer,ifirst,ilength)
int *unit;
unsigned long int *ifirst,*ilength;
char *buffer;
{
  unsigned long int length;
  PSI_FPTR first;

  first = (PSI_FPTR) (*ifirst-1)*sizeof(int);
  length = (*ilength)*sizeof(int);

  current_unit = *unit;

  if (first + length < first) {
    fprintf(stderr,"iordr_: *unit=%d, first=%lu, length=%lu\n",
            *unit, first, length);
    io_q_abort();
    }

  if (ud[*unit].method == PSIIO_SEQUENTIAL) {
    sequential_iordr(ud[*unit].ptr.sequential,buffer,first,length);
    }
  else if (ud[*unit].method == PSIIO_R_ASYNC) {
    r_async_iordr(ud[*unit].ptr.r_async,buffer,first,length);
    }
  else if (ud[*unit].method == PSIIO_S_ASYNC) {
    s_async_iordr(ud[*unit].ptr.s_async,buffer,first,length);
    }
  else if (ud[*unit].method == PSIIO_RAM) {
    ram_iordr(ud[*unit].ptr.ram,buffer,first,length);
    }
  else if (ud[*unit].method == PSIIO_FORTRAN) {
    C_IOFR (unit,buffer,&first,&length);
    }
  else {
    fprintf(stderr,"iordr_: invalid method (%d) for unit %d\n",
            ud[*unit].method,*unit);
    io_q_abort();
    }
  }

void
ioabort()
{
  fprintf(stderr,"ioabort: current unit = %d\n",current_unit);
  }

/* call ioabort to print out the unit number and then call the
 * fortran abort routine to flush out all of the buffers. */
void
io_q_abort()
{
  ioabort();

  /* Now call the fortran abort routine. */
#if defined(AIX) || defined(APOLLO)
  qabort();
#else
  qabort_();
#endif

  }
