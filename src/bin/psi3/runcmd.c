#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>

void runcmd(int *errcod, char *cmd)
{

  /*
  i = system(cmd);
  *errcod = WEXITSTATUS(i);
  */

  /* this version has the advantage that you can ctrl-C it */
  if (fork() == 0) {
    execl("/bin/sh","/bin/sh","-c",cmd,NULL);
    fprintf(stderr,"psi: execl returned?? Is the /bin/sh command missing?\n");
    exit(1);
  }

  wait(errcod);

  /*
  fprintf(stderr, "errcod in runcmd before filter is %d\n", i);
  *errcod = WEXITSTATUS(i);
  */
}

