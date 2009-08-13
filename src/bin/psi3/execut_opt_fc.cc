/*! \file
    \ingroup PSI3
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <psifiles.h>

extern "C" FILE *outfile;

namespace psi { namespace psi3 {

extern int execut(char **module_names, int num_modules, int depth);
extern void psi3_abort(void);
extern void runcmd(int *errcod, char *cmd);

int execut_opt_fc(char **exec, int nexec, char **exec_fc, int nexec_fc)
{
  int i, j, nrep, inloop, nopt, I;
  bool balked_last_time = false;
  int errcod, start_loop;

  // compute force constants
  execut(exec_fc,nexec_fc,0);

  // do geometry optimization prep
  for (i=0; i<nexec; i++) {

    if (!strcmp(exec[i], "REPEAT")) {
       sscanf(exec[++i],"%d",&nopt);
       start_loop = ++i;
       break;
    }

    fprintf(stderr," %s\n",exec[i]);
    runcmd(&errcod,exec[i]);
    //fprintf(stderr," errcod after runcmd is %d\n", errcod);

    if (WIFSIGNALED(errcod)) {
      fprintf(stdout,"\nCommand %s was terminated with signal %d", exec[i], WTERMSIG(errcod));
      psi3_abort();
    }
    errcod = WEXITSTATUS(errcod);
    //fprintf(stderr,"  errcod after wexitstatus is %d\n",errcod);
  }

  // do geometry optimization loop
  for (I=0 ; I<nopt; ++I) {
    fprintf(stderr," CYCLE %d of %d\n", I+1, nopt);

    for (j=start_loop; j<nexec; ++j) {

      fprintf(stderr," %s\n",exec[j]);
      runcmd(&errcod,exec[j]);
      //fprintf(stderr," errcod after runcmd is %d\n", errcod);

      if (WIFSIGNALED(errcod)) {
        fprintf(stdout,"\nCommand %s was terminated with signal %d", exec[j], WTERMSIG(errcod));
        psi3_abort();
      }
      errcod = WEXITSTATUS(errcod);

      if (!strcmp(exec[j], "optking --opt_step")) {

    //fprintf(stderr,"  errcod after wexitstatus is %d\n",errcod);

        if (errcod == PSI_RETURN_FAILURE) { //failed
          psi3_abort();
        }
        else if (errcod == PSI_RETURN_BALK) { //balked
          if (balked_last_time) {
            fprintf(outfile,"\nBalked two times consecutively. Aborting.\n");
            fprintf(outfile,"Unable to take an acceptable step after computing force constants.\n");
            fprintf(stdout,"\nBalked two times consecutively. Aborting.\n");
            fprintf(stdout,"Unable to take an acceptable step after computing force constants.\n");
            psi3_abort();
          }
          else balked_last_time = true;

          execut(exec_fc,nexec_fc,1);
          I--;
        }
        else if (errcod == PSI_RETURN_ENDLOOP) { //converged
          I=nopt;
        }
        else { // took step
          balked_last_time = false;
        }

        break; // skip 'end' and 'psiclean'
      }
    }
  }

  return(i);
}

}} // namespace psi::psi3
