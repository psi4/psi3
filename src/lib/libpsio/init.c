/*!
   \file init.c
   \ingroup (PSIO)
*/

#include <stdio.h>
#include <stdlib.h>
#include <libpsio/psio.h>

/* Definition of global data */
psio_lib* _default_psio_lib_ = NULL;
int _psio_error_exit_code_ = 1;
psio_address PSIO_ZERO = {0,0};

extern int psio_ipv1_config();

int __psio_init(psio_lib* Lib)
{
  int i,j;
  char *userhome;
  FILE *psirc;

  Lib->psio_unit = (psio_ud *) malloc(sizeof(psio_ud)*PSIO_MAXUNIT);
#ifdef PSIO_STATS
  Lib->psio_readlen = (ULI *) malloc(sizeof(ULI) * PSIO_MAXUNIT);
  Lib->psio_writlen = (ULI *) malloc(sizeof(ULI) * PSIO_MAXUNIT);
#endif
  Lib->state = 1;

  /* if in PSI3 -- configure libpsio using libipv1 */
  psio_ipv1_config();

  if(Lib->psio_unit == NULL) {
    fprintf(stderr, "Error in PSIO_INIT()!\n");
    exit(_psio_error_exit_code_);
  }

  for(i=0; i < PSIO_MAXUNIT; i++) {
#ifdef PSIO_STATS
    Lib->psio_readlen[i] = Lib->psio_writlen[i] = 0;
#endif      
    Lib->psio_unit[i].numvols = 0;
    for(j=0; j < PSIO_MAXVOL; j++) {
      Lib->psio_unit[i].vol[j].path = NULL;
      Lib->psio_unit[i].vol[j].stream = -1;
    }
    Lib->psio_unit[i].toclen = 0;
    Lib->psio_unit[i].toc = NULL;
  }

  return(1);
}

/*!
** PSIO_INIT(): Allocates global memory needed by the I/O routines.
**
** No arguments.
**
** \ingroup (PSIO)
*/

int psio_init(void)
{
  if (!_default_psio_lib_) {
    _default_psio_lib_ = (psio_lib *)malloc(sizeof(psio_lib));
    if (_default_psio_lib_ == NULL) {
      fprintf(stderr,"LIBPSIO::init() -- failed to allocate the memory");
      exit(_psio_error_exit_code_);
    }
  }

  return __psio_init(_default_psio_lib_);
}

/*!
** PSIO_STATE(): Returns state of the library (1=initialized, 0=noninitialized).
**
** No arguments.
**
** \ingroup (PSIO)
*/

int psio_state()
{
  return _default_psio_lib_->state;
}

