
#include <psiconfig.h>

#include <stdio.h>
#include <string.h>
#include <time.h>

#ifdef HAVE_SYS_TIMES_H
#include <sys/times.h>
#endif

#ifdef HAVE_SYS_PARAM_H
#include <sys/param.h>
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifndef HZ
#define HZ 60
#endif

#if FCLINK==1
#define C_TIMNOW timnow_
#define C_CFDATE cfdate_
#define C_GTHOST gthost_
#define C_CTIMES ctimes_
#elif FCLINK==2 
#define C_TIMNOW timnow
#define C_CFDATE cfdate
#define C_GTHOST gthost
#define C_CTIMES ctimes
#else
#define C_TIMNOW TIMNOW
#define C_CFDATE CFDATE
#define C_GTHOST GTHOST
#define C_CTIMES CTIMES
#endif

void
C_TIMNOW(psi_int *now)
{
  *now = (psi_int) time(NULL);
}

void
C_CFDATE(psi_int *now, char fdate[80])
{
  int i;

  strncpy(fdate,ctime((time_t*)now),80);

  /* get rid of newline */
  fdate[strlen(fdate)-1] = '\0';

  /* and pad for fortran */
  for (i=strlen(fdate); i < 80; i++) {
    fdate[i] = ' ';
  }
}

void
C_GTHOST(char hostnm[80])
{
  int i;

  if (gethostname(hostnm,80) == -1) {
    strcpy(hostnm,"unknown");
  }

  /* pad hostnm */
  for (i=strlen(hostnm); i < 80; i++) {
    hostnm[i] = ' ';
  }
}

void
C_CTIMES(psi_real *utime, psi_real *stime)
{
  struct tms t;

  times(&t);

  *utime = (double)t.tms_utime/(double)HZ;
  *stime = (double)t.tms_stime/(double)HZ;
}

