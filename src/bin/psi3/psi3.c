/*
 * PSI3
 * 
 * New PSI driver
 *
 * This is the program that executes all the PSI modules.  It decides
 * what modules to run by parsing the file [psi-binary-path]/share/psi.dat.
 * 
 * C. David Sherrill
 * Georgia Institute of Technology
 * July 2002
 */

#include <stdio.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <psifiles.h>          /* where return values are */
#include <sys/types.h>
#include <sys/wait.h>

#define MXEXEC 100
#define MAX_EXEC_STR 80

FILE *infile, *outfile;

void psi3_abort(void);
int execut(char **module_names, int num_modules, int depth);
extern char **parse_var(int *nvars, int mxvars, char *name);
extern void runcmd(int *errcod, char *cmd);
void parse_cmdline(int argc, char *argv[]);
  

int main(int argc, char *argv[])
{
  FILE *psidat;
  char *wfn, *dertyp, *reftyp, *calctyp, **exec, proced[132];
  int check;
  int i,j,nexec=0,rdepth=0;
  int errcod;

  enum CalcCode {
    SP,               /* single-point */
    OPT,              /* optimization */
    DISP,             /* displacements */
    FREQ,             /* frequencies */
    SYMM_FREQ         /* frequencies for symmetric modes */
  } JobType;

  parse_cmdline(argc,argv);

  fprintf(outfile, "\nThe PSI3 Execution Driver\n");

  psidat = fopen(SHARE, "r");

  ip_set_uppercase(1);
  ip_initialize(infile, outfile);
  ip_append(psidat,outfile);
  ip_cwk_clear();
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(":PSI");
  fclose(psidat);


  /* read in input parameters */
  errcod = ip_string("WFN",&wfn,0);
  if (errcod != IPE_OK) {
    fprintf(outfile,"Error: Could not read required keyword WFN from input\n");
    psi3_abort();
  }

  errcod = ip_string("REFERENCE",&reftyp,0);
  if (errcod == IPE_KEY_NOT_FOUND) {
    reftyp = (char *) malloc(sizeof(char)*4);
    strcpy(reftyp,"RHF");
  }

  errcod = ip_boolean("CHECK",&check,0);
  if (errcod == IPE_KEY_NOT_FOUND)
    check = 0;
  else if (errcod != IPE_OK) {
    fprintf(outfile, 
            "Error: A problem arose reading the optional keyword 'check'.\n");
    psi3_abort();
  }

  
  /* this is the old way to get the calculation type
  errcod = ip_boolean("DISP",&disp,0);
  if (errcod == IPE_KEY_NOT_FOUND)
    disp = 0;
  else if (errcod != IPE_OK) {
    fprintf(outfile, 
            "Error: A problem arose reading the optional keyword 'disp'.\n");
    psi3_abort();
  }

  errcod = ip_boolean("OPT",&opt,0);
  if (errcod == IPE_KEY_NOT_FOUND)
    opt = 0;
  else if (errcod != IPE_OK) {
    fprintf(outfile, 
            "Error: A problem arose reading the optional keyword 'opt'.\n");
    psi3_abort();
  }

  errcod = ip_boolean("FREQ",&freq,0);
  if (errcod == IPE_KEY_NOT_FOUND)
    freq = 0;
  else if (errcod != IPE_OK) {
    fprintf(outfile, 
            "Error: A problem arose reading the optional keyword 'freq'.\n");
    psi3_abort();
  }
  */

  /* get the calculation type */
  errcod = ip_string("JOBTYPE",&calctyp,0);
  if (errcod == IPE_KEY_NOT_FOUND) {
    calctyp = (char *) malloc(sizeof(char)*3);
    strcpy(calctyp,"SP");
  }
  else if (errcod != IPE_OK) {
    fprintf(outfile, "Trouble reading JOBTYPE from input file\n");
    psi3_abort();
  }

  /* set JobType depending on JOBTYPE */
  if ((strcmp(calctyp,"FREQ")==0) || (strcmp(calctyp,"FREQUENCY")==0))
    JobType = FREQ;
  else if ((strcmp(calctyp,"OPT")==0) || (strcmp(calctyp,"OPTIMIZATION")==0))
    JobType = OPT;
  else if ((strcmp(calctyp,"DISP")==0) || (strcmp(calctyp,"DISPLACEMENTS")==0))
    JobType = DISP;
  else if ((strcmp(calctyp,"FREQ_SYMM")==0)||(strcmp(calctyp,"FREQ-SYMM")==0)||
           (strcmp(calctyp,"SYMM_FREQ")==0)||(strcmp(calctyp,"SYMM-FREQ")==0)||
	   (strcmp(calctyp,"SYMMETRIC_FREQUENCY")==0)) JobType = SYMM_FREQ;
  else if ((strcmp(calctyp,"SP")==0) || (strcmp(calctyp,"SINGLE-POINT")==0) ||
           (strcmp(calctyp,"SINGLE_POINT")==0) ||
	   (strcmp(calctyp,"FORCE")==0)) JobType = SP;
  else {
    fprintf(outfile,"Error: Unrecognized calculation type %s\n", calctyp);
    psi3_abort();
  }

  /* default dertyp depends on calculation type */
  errcod = ip_string("DERTYPE",&dertyp,0);
  if (errcod == IPE_KEY_NOT_FOUND) {
    dertyp = (char *) malloc(sizeof(char)*10);
    if (strcmp(calctyp,"FORCE")==0) 
      strcpy(dertyp,"FIRST");
    else if (JobType == FREQ || JobType == SYMM_FREQ) {
      if (strcmp(wfn,"SCF")==0) strcpy(dertyp,"SECOND");
      else strcpy(dertyp,"FIRST"); /* guess analyt grads unless overridden */
    }
    else
      strcpy(dertyp,"NONE");
  }


  /* make some basic checks on the requested computation type */
  if ((strcmp(reftyp,"RHF")!=0) && (strcmp(reftyp,"ROHF")!=0) &&
      (strcmp(reftyp,"UHF")!=0) && (strcmp(reftyp,"TWOCON")!=0)) 
  {
    fprintf(outfile,"Error: bad 'reference'.\n");
    fprintf(outfile,"Must be one of: RHF, ROHF, UHF, TWOCON\n");
    psi3_abort();
  }

  if ((strcmp(dertyp,"NONE")!=0) && (strcmp(dertyp,"FIRST")!=0) &&
      (strcmp(dertyp,"SECOND")!=0))
  {
    fprintf(outfile,"Error: bad 'dertype'.\n");
    fprintf(outfile,"Must be one of: NONE, FIRST, SECOND\n");
    psi3_abort();
  }

  /* print out what type of calculation it is */
  fprintf(outfile, "\nPSI3 will perform a %s %s ", reftyp, wfn);
  if (JobType == OPT) {
    if (strcmp(dertyp,"FIRST")==0) 
      fprintf(outfile,"optimization via analytic gradients.\n");
    else if (strcmp(dertyp,"NONE")==0)
      fprintf(outfile,"optimization via energy points.\n");
    else 
      fprintf(outfile,"optimization using DERTYPE=%s.\n", dertyp);
  }
  else if (JobType == DISP) {
    if (strcmp(dertyp,"NONE")==0)
      fprintf(outfile,"displacement computation.\n");
    else if (strcmp(dertyp,"FIRST")==0)
      fprintf(outfile, "gradient displacement computation.\n");
    else if (strcmp(dertyp,"SECOND")==0)
      fprintf(outfile, "Hessian displacement computation.\n");
  }
  else if (JobType == FREQ) {
    if (strcmp(dertyp,"NONE")==0)
      fprintf(outfile,"frequency computation via energies.\n");
    else if (strcmp(dertyp,"FIRST")==0)
      fprintf(outfile,"frequency computation via gradients.\n");
    else if (strcmp(dertyp,"SECOND")==0)
      fprintf(outfile, "analytic frequency computation.\n");
  }
  else if (JobType == SYMM_FREQ) {
    if (strcmp(dertyp,"NONE")==0)
      fprintf(outfile,"symmetric frequency computation via energies.\n");
    else if (strcmp(dertyp,"FIRST")==0)
      fprintf(outfile,"symmetric frequency computation via gradients.\n");
    else if (strcmp(dertyp,"SECOND")==0)
      fprintf(outfile, "analytic symmetric frequency computation.\n");
  }
  else if (JobType == SP) { 
    if (strcmp(dertyp,"NONE")==0) fprintf(outfile, "energy");
    else if (strcmp(dertyp,"FIRST")==0) fprintf(outfile, "gradient");
    else if (strcmp(dertyp,"SECOND")==0) fprintf(outfile, "Hessian");
    else fprintf(outfile, "unrecognized-dertype");
    fprintf(outfile, " computation.\n");
  }
  else { 
    fprintf(outfile, "calculation of unrecognized type.\n");
  }

  if (ip_exist("EXEC",0)) { /* User-specified EXEC statement */
    fprintf(outfile,"Using the user provided execution list from 'exec'.\n"); 
    exec = parse_var(&nexec,MXEXEC,"EXEC");
  }
  else { /* Default sequence of modules from psi.dat */
    /* construct the name of the procedure from dertyp, reftyp, and wfn */
    strcpy(proced,"");
    strcat(proced,wfn);
    strcat(proced,reftyp);
    if (strcmp(dertyp,"NONE")==0) strcat(proced,"ENERGY");
    else strcat(proced,dertyp);
    if (JobType==OPT) strcat(proced,"OPT");
    else if (JobType==FREQ) strcat(proced,"FREQ");
    else if (JobType==SYMM_FREQ) strcat(proced,"SYMM_FREQ");
    else if (JobType==DISP) strcat(proced,"DISP");

    fprintf(outfile, "Calculation type string = %s\n", proced);

    if (!ip_exist(proced,0)) {
      fprintf(outfile,"Error: Did not find a valid calculation type,\n");
      fprintf(outfile,
              "probably because of incompatible wfn,dertype,reference.\n");
      fprintf(outfile,"Full calculation type is: %s\n", proced); 
      psi3_abort();
    }
    exec = parse_var(&nexec,MXEXEC,proced);
    /* fprintf(outfile, "nexec = %d\n", nexec); */
  }

  if (check) {
    fprintf(outfile,"\n'CHECK' is YES, so nothing will be executed.\n");
    fprintf(outfile,"The following programs would otherwise be executed:\n");
  }
  else
    fprintf(outfile,"The following programs will be executed:\n");

  rdepth = 0;
  for (i=0; i<nexec; i++) {
    if (strcmp(exec[i],"END")==0) rdepth=rdepth-1;
    for (j=0; j<rdepth; j++) fprintf(outfile," ");
    if (strcmp(exec[i],"REPEAT")==0) {
      rdepth=rdepth+1;
      fprintf(outfile," %s %s\n",exec[i],exec[i+1]);
      i++;
    }
    else 
      fprintf(outfile," %s\n",exec[i]);
  }

  fprintf(outfile,"\n");
  if (!check) execut(exec,nexec,0); 


  /* clean up and free memory */
  ip_done();
  free(wfn);
  free(dertyp);
  free(reftyp);
  for (i=0; i<MXEXEC; i++) free(exec[i]);

}


void psi3_abort(void)
{
  fprintf(outfile,"\nPSI3 exiting.\n");
  exit(0);
}

int execut(char **exec, int nexec, int depth)
{
  int i, j, nrep, inloop;
  int errcod;
  char spaces[80];

  sprintf(spaces, "");
  for (i=0; i<depth; i++)
    strcat(spaces, " ");

  for (i=0; i<nexec; i++) {
    if (strcmp(exec[i],"END")==0) return(i);
    if (strcmp(exec[i],"REPEAT")!=0) {
      fprintf(outfile,"%s%s\n", spaces, exec[i]);
      runcmd(&errcod,exec[i]);

      /* fprintf(stderr,"%serrcod before filter is %d\n",spaces,errcod); */

      /* exited with signal */
      if (WIFSIGNALED(errcod)) {
        fprintf(outfile,"\nCommand was terminated with signal %d",
                WTERMSIG(errcod));
	psi3_abort();
      }

      errcod = WEXITSTATUS(errcod);

      /* fprintf(stderr,"%serrcod after filter is %d\n",spaces,errcod); */

      /* if we're getting ndisp from optking */
      if (strcmp(exec[i],"optking --disp_symm")==0 ||
          strcmp(exec[i],"optking --disp_all")==0) {
        /* fprintf(outfile, "optking got exit code %d\n", errcod); */
        for (j=i; j<nexec; j++) {
          /* fprintf(outfile,"Scanning exec[%d] = %s\n", j, exec[j]);  */
          if (strcmp(exec[j],"NUM_DISP")==0) {	
            sprintf(exec[j],"%d",errcod);
            break;
          }
        }
	/* if we didn't find it, assume was already replaced in prev loop */
      }

      /* if we need to break out of the loop we are in */
      else if (errcod == PSI_RETURN_ENDLOOP) {
        for (j=i; j<nexec; j++) {
          if (strcmp(exec[j],"END")==0) break;
	}
	if (j==nexec) {
          fprintf(outfile,"PSI3 caught ENDLOOP signal but no END below!\n");
          psi3_abort();
        }
        else return(-j);
      }
      /* if we got a 'command failed' failed flag */
      else if (errcod == PSI_RETURN_FAILURE) {
        fprintf(outfile,"\nCommand has returned a fail status.");
	psi3_abort();
      }
    }
    else {
      sscanf(exec[i+1], "%d", &nrep);
      fprintf(outfile, "%sREPEAT (%d)\n", spaces, nrep);
      for (j=0; j<nrep; j++) {
        fprintf(outfile, "%s CYCLE %d\n", spaces, j+1);
        inloop = execut(&(exec[i+2]),nexec-i-2,depth+1);
	if (inloop < 0) {
          inloop = -inloop;
          break;
        }
      }
      i = i + inloop + 2;
    }
  }
  return(i);
}

/* this code is essentially the same as in the psi_start() lib funct */
int parse_cmdline(int argc, char *argv[])
{
  int i;
  int found_if_np = 0;          /* found input file name without -i */
  int found_of_np = 0;          /* found output file name without -o */
  int found_fp_np = 0;          /* found file prefix name without -p */
  int found_if_p = 0;           /* found input file name with -i */
  int found_of_p = 0;           /* found output file name with -o */
  int found_fp_p = 0;           /* found file prefix name with -p */
  char *ifname, *ofname, *fprefix, *arg;

  /* process command-line arguments in sequence */
  for(i=0; i<argc; i++) {
    arg = argv[i];
    if (!strcmp(arg,"-f") && !found_if_p) {
      ifname = argv[++i];
      found_if_p = 1;
    }
    else if (!strcmp(arg,"-o") && !found_of_p) {
      ofname = argv[++i];
      found_of_p = 1;
    }
    else if (!strcmp(arg,"-p") && !found_fp_p) {
      fprefix = argv[++i];
      found_fp_p = 1;
    }
    else if (arg[0] == '-') {
      fprintf(stderr, "Error: unrecognized command-line argument %s\n", arg);
      return(0);
    }
    else if (!found_if_np) {
      ifname = arg;
      found_if_np = 1;
    }
    else if (!found_of_np) {
      ofname = arg;
      found_of_np = 1;
    }
    else if (!found_fp_np) {
      fprefix = arg;
      found_fp_np = 1;
    }
    else {
      fprintf(stderr, "Error: too many command-line arguments given\n");
      return(0);
    }
  }
  

  /* check if some arguments were specified in both prefixed and nonprefixed 
   * form */
  if (found_if_p && found_if_np) {
    fprintf(stderr, 
            "Error: input file name specified both with and without -f\n");
    fprintf(stderr, 
            "Usage: (module) [options] -f input -o output [-p prefix]  OR\n");
    fprintf(stderr, "       (module) [options] input output [prefix]\n");
    return(0);
  }
  if (found_of_p && found_of_np) {
    fprintf(stderr, 
            "Error: output file name specified both with and without -o\n");
    fprintf(stderr, 
            "Usage: (module) [options] -f input -o output [-p prefix]  OR\n");
    fprintf(stderr, "       (module) [options] input output [prefix]\n");
    return(0);
  }
  if (found_fp_p && found_fp_np) {
    fprintf(stderr, 
            "Error: file prefix specified both with and without -p\n");
    fprintf(stderr, 
            "Usage: (module) [options] -f input -o output -p prefix  OR\n");
    fprintf(stderr, "       (module) [options] input output prefix\n");
    return(0);
  }

 
  /* set the environmental variables the modules will look for */ 
  if (ifname != NULL)
    setenv("PSI_INPUT",ifname);
  if (ofname != NULL)
    setenv("PSI_OUTPUT",ofname);
  if (fpname != NULL)
    setenv("PSI_PREFIX",fpname);

}

