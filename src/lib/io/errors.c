
#include "includes.h"
#include "param.h"
#include "types.h"
#include "errno.h"

void
no_path_given(name)
char *name;
{
  fprintf(stderr,"%s: no path given\n",name);
  io_q_abort();
  }

void
malloc_check(caller,data)
char *caller;
char *data;
{
  if (!data) {
    fprintf(stderr,"%s: malloc failed\n");
#if defined(AIX)
    fprintf(stderr,"errno = %d\n",errno);
#else
    perror("malloc");
#endif
    io_q_abort();
    }
  }

void
fopen_check(caller,path,data)
char *caller;
char *path;
char *data;
{
  if (!data) {
    fprintf(stderr,"%s: fopen failed for %s\n",caller,path);
#if defined(AIX)
    fprintf(stderr,"errno = %d\n",errno);
#else
    perror("fopen");
#endif
    io_q_abort();
    }
  }

void
fread_error(caller)
char *caller;
{
  fprintf(stderr,"%s: fread failed\n",caller);
#if defined(AIX)
  fprintf(stderr,"errno = %d\n",errno);
#else
  perror("fread");
#endif
  io_q_abort();
  }

void
fwrite_error(caller)
char *caller;
{
  fprintf(stderr,"%s: fwrite failed\n",caller);
#if defined(AIX)
  fprintf(stderr,"errno = %d\n",errno);
#else
  perror("fwrite");
#endif
  io_q_abort();
  }
