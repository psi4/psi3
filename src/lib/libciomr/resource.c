#include "includes.h"
#include "common.h"

#define STRING_LENGTH 80

#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#if defined(sun)
# include <unistd.h>
#endif

void
resource_command()
{
  char string[STRING_LENGTH];
  struct rusage r;
  FILE *output;

  output=outfile;

  getrusage(RUSAGE_SELF,&r);

      print_time('u',r.ru_utime.tv_sec,r.ru_utime.tv_usec);
      print_time('s',r.ru_stime.tv_sec,r.ru_stime.tv_usec);
      fprintf(output,"  maxrss: %d kB\n",(r.ru_maxrss*getpagesize())/1000);
      fprintf(output,"  ixrss: %d\n",r.ru_ixrss);
      fprintf(output,"  idrss: %d\n",r.ru_idrss);
      fprintf(output,"  isrss: %d\n",r.ru_isrss);
      fprintf(output,"  minflt: %d\n",r.ru_minflt);
      fprintf(output,"  majflt: %d\n",r.ru_majflt);
      fprintf(output,"  nswap: %d\n",r.ru_nswap);
      fprintf(output,"  inblock: %d\n",r.ru_inblock);
      fprintf(output,"  oublock: %d\n",r.ru_oublock);
      fprintf(output,"  nvcsw: %d\n",r.ru_nvcsw);
      fprintf(output,"  nivcsw: %d\n",r.ru_nivcsw);

  fflush(outfile);

  }

print_time(ch,sec,usec)
   char ch;
   long sec, usec;
{
  int i;
  char usec_string[7];

  sprintf(usec_string,"%6ld",usec);
  for (i=0; i<6; i++) {
    if (usec_string[i] == ' ') usec_string[i] = '0';
    }

  fprintf(outfile,"  %ctime: %ld.%s\n",
           ch,sec,usec_string);
  }

