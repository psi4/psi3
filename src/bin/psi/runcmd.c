
#include <stdio.h>

void
#if defined(APOLLO) || defined(AIX)
runcmd(errcod,cmd,length)
#else
runcmd_(errcod,cmd,length)
#endif
int *errcod;
char *cmd;
int *length;
{
  int i,last;
  char copy[200];

  if (*length >= 200) {
    fprintf(stderr,"INTERNAL ERROR: *length >= 200");
    exit(1);
    }

  last = -1;
  for (i=0; i<*length; i++) {
    if (cmd[i]!=' ') last = i;
    copy[i] = cmd[i];
    }
  copy[last+1] = '\0';

  /* Due to AIX's bug in the system lib routine this won't work: */
  /* *errcod = system(copy); */
  if (fork() == 0) {
    execl("/bin/sh","/bin/sh","-c",copy,NULL);
    fprintf(stderr,"psi: execl returned?? Is the /bin/sh command missing?\n");
    exit(1);
    }
  /* Wait for the child to die. */
  wait(errcod);
  }

