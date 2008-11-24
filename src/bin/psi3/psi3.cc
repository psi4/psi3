/*! \file
    \ingroup PSI3
    \brief This is the main PSI3 driver program
*/
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

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <cstring>
#include <ctime>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <psifiles.h>          /* where return values are */
#include <sys/types.h>
#include <sys/wait.h>
#include <libpsio/psio.h>
#include <psifiles.h>

#include <sys/times.h>

#define MXEXEC 100
#define MAX_EXEC_STR 80

#if defined HAVE_DECL_SETENV && !HAVE_DECL_SETENV
  extern int setenv(const char *, const char *, int);
#endif

extern "C" {
FILE *infile, *outfile;
char *psi_file_prefix;
}

namespace psi { namespace psi3 {

time_t time_start, time_end;
struct tms total_tmstime;


int get_ndisp(void);
/* 
  the following must stay in scope throughout the run, else the 
  environmental variables will vanish!
*/

char tmpstr_input[200], tmpstr_output[200], tmpstr_prefix[200];
const char *ofname = NULL;

void psi3_abort(void);
int execut(char **module_names, int num_modules, int depth);
extern char **parse_var(int *nvars, int mxvars, const char *name);
extern void runcmd(int *errcod, char *cmd);
int parse_cmdline(int argc, char *argv[]);
void print_welcome(FILE *outfile);
void tstart_psi3driver(FILE *outfile);
void tstop_psi3driver(FILE *outfile);

/* boolean for automatically running input program first in all procedures */
int auto_input;
/* boolean for checking input w/o running any programs */
int auto_check;
/* boolean for whether called from dboc module (if yes - ignore jobtype=dboc) */
int called_from_dboc;
/* boolean for whether $done should be called */
int call_done;
/* boolean to know if the driver should write to the output file or not 
   (the psi test driver wants to get stdout info only by running
   auto_check, and in such a case we don't want to write to output...) */
int have_outfile;
}} // namespace psi::psi3

int main(int argc, char *argv[])
{
  using namespace psi::psi3;
  FILE *psidat;
  char *psidat_dirname, *psidat_filename;
  char *wfn, *dertyp, *reftyp, *calctyp, *jobtype, **exec, proced[132];
  char tmpstr[133];
  int check=0;
  int i,j,nexec=0,rdepth=0;
  int direct=0;
  int errcod;
  int plus_d=0;
  char **input_exec;
  int nexec_input;

  enum CalcCode {
    SP,           /* single-point */
    OPT,          /* optimization */
    DISP,         /* displacements */
    FREQ,         /* frequencies */
    SYMM_FC,      /* force constants in internal coordinates, symmetric modes */
    FC,           /* force constants in internal coordinates, all modes */
    OEPROP,       /* one-electron properties */
    DBOC          /* compute Diagonal Born-Oppenheimer Correction (DBOC) 
                     by finite difference */
  } JobType;

  // Set these to a known value
  infile = NULL;
  outfile = NULL;

  if (!parse_cmdline(argc,argv))
    psi3_abort();

  if (have_outfile) {
    tstart_psi3driver(outfile);
    print_welcome(outfile);
  }

  fprintf(stdout, "\n\n The PSI3 Execution Driver \n");

  /* To find psi.dat first check the environment, then its location 
     after installation */
  psidat_dirname = getenv("PSIDATADIR");
  if (psidat_dirname != NULL) {
    char* tmpstr = (char *) malloc(sizeof(char)*(strlen(psidat_dirname)+9));
    sprintf(tmpstr,"%s/psi.dat",psidat_dirname);
    psidat_filename = tmpstr;
  }
  else
    psidat_filename = strdup(INSTALLEDPSIDATADIR "/psi.dat");
  psidat = fopen( psidat_filename, "r");
  if (psidat == NULL) {
    fprintf(stdout, "Warning: didn't find psi.dat at %s\n", psidat_filename);
  }
  free(psidat_filename);

  ip_set_uppercase(1);
  ip_initialize(infile, stdout);
  ip_append(psidat,stdout);
  ip_cwk_clear();
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(":PSI");
  fclose(psidat);


  /* read in input parameters */
  errcod = ip_string("WFN",&wfn,0);
  if (errcod != IPE_OK) {
    fprintf(stdout,"Error: Could not read required keyword WFN from input\n");
    psi3_abort();
  }

  errcod = ip_boolean("EMPIRICAL_DISPERSION",&plus_d,0);
  if (plus_d) {
    if (strcmp(wfn,"SCF")==0) {
      free(wfn);
      wfn = strdup("SCF+D");  
    }
    else {
      fprintf(outfile, "Error: EMPIRICAL_DISPERSION not supported for\n");
      fprintf(outfile, " wavefunctions other than wfn=SCF.\n");
    }
  }

  errcod = ip_boolean("DIRECT",&direct,0);

    
  errcod = ip_string("REFERENCE",&reftyp,0);
  if (errcod == IPE_KEY_NOT_FOUND) {
    reftyp = (char *) malloc(sizeof(char)*4);
    strcpy(reftyp,"RHF");
  }

  if(!auto_check) {
    errcod = ip_boolean("CHECK",&check,0);
    if (errcod == IPE_KEY_NOT_FOUND)
      check = 0;
    else if (errcod != IPE_OK) {
      fprintf(stdout, 
              "Error: A problem arose reading the optional keyword 'check'.\n");
      psi3_abort();
    }
  }

  
  /* this is the old way to get the calculation type
     errcod = ip_boolean("DISP",&disp,0);
     if (errcod == IPE_KEY_NOT_FOUND)
     disp = 0;
     else if (errcod != IPE_OK) {
     fprintf(stdout, 
     "Error: A problem arose reading the optional keyword 'disp'.\n");
     psi3_abort();
     }

     errcod = ip_boolean("OPT",&opt,0);
     if (errcod == IPE_KEY_NOT_FOUND)
     opt = 0;
     else if (errcod != IPE_OK) {
     fprintf(stdout, 
     "Error: A problem arose reading the optional keyword 'opt'.\n");
     psi3_abort();
     }

     errcod = ip_boolean("FREQ",&freq,0);
     if (errcod == IPE_KEY_NOT_FOUND)
     freq = 0;
     else if (errcod != IPE_OK) {
     fprintf(stdout, 
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
    fprintf(stdout, "Trouble reading JOBTYPE from input file\n");
    psi3_abort();
  }

  /* set JobType depending on JOBTYPE */
  if ((strcmp(calctyp,"FREQ")==0) || (strcmp(calctyp,"FREQUENCY")==0))
    JobType = FREQ;
  else if ((strcmp(calctyp,"OPT")==0) || (strcmp(calctyp,"OPTIMIZATION")==0))
    JobType = OPT;
  else if ((strcmp(calctyp,"DISP")==0) || (strcmp(calctyp,"DISPLACEMENTS")==0))
    JobType = DISP;
  else if ( (strcmp(calctyp,"SYMM_FC")==0)||(strcmp(calctyp,"SYMM-FC")==0))
    JobType = SYMM_FC;
  else if ( (strcmp(calctyp,"FC")==0) )
    JobType = FC;
  else if ((strcmp(calctyp,"SP")==0) || (strcmp(calctyp,"SINGLE-POINT")==0) ||
           (strcmp(calctyp,"SINGLE_POINT")==0) ||
	   (strcmp(calctyp,"FORCE")==0)) JobType = SP;
  else if ((strcmp(calctyp,"OEPROP")==0))
    JobType = OEPROP;
  else if ((strcmp(calctyp,"DBOC")==0) || strcmp(calctyp,"BODC")==0) {
    if (called_from_dboc)
      JobType = SP;
    else
      JobType = DBOC;
  }
  else {
    fprintf(stdout,"Error: Unrecognized calculation type %s\n", calctyp);
    psi3_abort();
  }

  /* default dertyp depends on calculation type */
  errcod = ip_string("DERTYPE",&dertyp,0);
  if (errcod == IPE_KEY_NOT_FOUND) {
    dertyp = (char *) malloc(sizeof(char)*10);
    if (strcmp(calctyp,"FORCE")==0) 
      strcpy(dertyp,"FIRST");
    else if (JobType == OPT) strcpy(dertyp, "FIRST");
    else if (JobType == FREQ || JobType == SYMM_FC) {
      if (strcmp(wfn,"SCF")==0) strcpy(dertyp,"SECOND");
      else strcpy(dertyp,"FIRST"); /* guess analyt grads unless overridden */
    }
    else
      strcpy(dertyp,"NONE");
  }
  /* DBOC doesn't need derivatives, only wave functions */
  if (JobType == DBOC || (JobType == SP && called_from_dboc)) {
    free(dertyp);
    dertyp = strdup("NONE");
  }

  /* make some basic checks on the requested computation type */
  if ((strcmp(reftyp,"RHF")!=0) && (strcmp(reftyp,"ROHF")!=0) &&
      (strcmp(reftyp,"UHF")!=0) && (strcmp(reftyp,"TWOCON")!=0))
    {
      fprintf(stdout,"Error: bad 'reference'.\n");
      fprintf(stdout,"Must be one of: RHF, ROHF, UHF, TWOCON\n");
      psi3_abort();
    }

  if ((strcmp(dertyp,"NONE")!=0) && (strcmp(dertyp,"FIRST")!=0) &&
      (strcmp(dertyp,"SECOND")!=0) && (strcmp(dertyp,"RESPONSE")!=0))
    {
      fprintf(stdout,"Error: bad 'dertype'.\n");
      fprintf(stdout,"Must be one of: NONE, FIRST, SECOND, or RESPONSE\n");
      psi3_abort();
    }

  /* print out what type of calculation it is */
  fprintf(stdout, "\nPSI3 will perform a %s %s ", reftyp, wfn);
  if (have_outfile)
  fprintf(outfile, "\nPSI3 will perform a %s %s ", reftyp, wfn);
  if (JobType == OPT) {
    if (strcmp(dertyp,"FIRST")==0) {
      fprintf(stdout,"optimization via analytic gradients.\n");
      if (have_outfile)
      fprintf(outfile,"optimization via analytic gradients.\n");
    }
    else if (strcmp(dertyp,"NONE")==0) {
      fprintf(stdout,"optimization via energy points.\n");
      if (have_outfile)
      fprintf(outfile,"optimization via energy points.\n");
    }
    else {
      fprintf(stdout,"optimization using DERTYPE=%s.\n", dertyp);
      if (have_outfile)
      fprintf(outfile,"optimization using DERTYPE=%s.\n", dertyp);
    }
  }
  else if (JobType == DISP) {
    if (strcmp(dertyp,"NONE")==0) {
      fprintf(stdout,"displacement computation.\n");
      if (have_outfile)
      fprintf(outfile,"displacement computation.\n");
    }
    else if (strcmp(dertyp,"FIRST")==0) {
      fprintf(stdout, "gradient displacement computation.\n");
      if (have_outfile)
      fprintf(outfile, "gradient displacement computation.\n");
    }
    else if (strcmp(dertyp,"SECOND")==0) {
      fprintf(stdout, "Hessian displacement computation.\n");
      if (have_outfile)
      fprintf(outfile, "Hessian displacement computation.\n");
    }
  }
  else if (JobType == FREQ) {
    if (strcmp(dertyp,"NONE")==0) {
      fprintf(stdout,"frequency computation via energies.\n");
      if (have_outfile)
      fprintf(outfile,"frequency computation via energies.\n");
    }
    else if (strcmp(dertyp,"FIRST")==0) {
      fprintf(stdout,"frequency computation via gradients.\n");
      if (have_outfile)
      fprintf(outfile,"frequency computation via gradients.\n");
    }
    else if (strcmp(dertyp,"SECOND")==0) {
      fprintf(stdout, "analytic frequency computation.\n");
      if (have_outfile)
      fprintf(outfile, "analytic frequency computation.\n");
    }
  }
  else if (JobType == SYMM_FC) {
    if (strcmp(dertyp,"NONE")==0) {
      fprintf(stdout,"symmetric frequency computation via energies.\n");
      if (have_outfile)
      fprintf(outfile,"symmetric frequency computation via energies.\n");
    }
    else if (strcmp(dertyp,"FIRST")==0) {
      fprintf(stdout,"symmetric frequency computation via gradients.\n");
      if (have_outfile)
      fprintf(outfile,"symmetric frequency computation via gradients.\n");
    }
    else if (strcmp(dertyp,"SECOND")==0) {
      fprintf(stdout, "analytic symmetric frequency computation.\n");
      if (have_outfile)
      fprintf(outfile, "analytic symmetric frequency computation.\n");
    }
  }
  else if (JobType == SP) { 
    if (strcmp(dertyp,"NONE")==0) {
      fprintf(stdout, "energy");
      if (have_outfile)
      fprintf(outfile, "energy");
    }
    else if (strcmp(dertyp,"FIRST")==0) {
      fprintf(stdout, "gradient");
      if (have_outfile)
      fprintf(outfile, "gradient");
    }
    else if (strcmp(dertyp,"SECOND")==0) {
      fprintf(stdout, "Hessian");
      if (have_outfile)
      fprintf(outfile, "Hessian");
    }
    else if (strcmp(dertyp,"RESPONSE")==0) {
      fprintf(stdout,"response property");
      if (have_outfile)
      fprintf(outfile,"response property");
    }
    else {
      fprintf(stdout, "unrecognized-dertype");
      if (have_outfile)
      fprintf(outfile, "unrecognized-dertype");
    }
    fprintf(stdout, " computation.\n");
    if (have_outfile)
    fprintf(outfile, " computation.\n");
  }
  else if (JobType == OEPROP) {
    fprintf(stdout, "one-electron properties computation.\n");
    if (have_outfile)
    fprintf(outfile, "one-electron properties computation.\n");
  }
  else if (JobType == DBOC) {
    fprintf(stdout,  
      "Diagonal Born-Oppenheimer Correction (DBOC) computation.\n");
    if (have_outfile)
    fprintf(outfile,  
      "Diagonal Born-Oppenheimer Correction (DBOC) computation.\n");
  }
  else { 
    fprintf(stdout, "calculation of unrecognized type.\n");
    if (have_outfile)
    fprintf(outfile, "calculation of unrecognized type.\n");
  }

  if (ip_exist("EXEC",0)) { /* User-specified EXEC statement */
    fprintf(stdout,"Using the user provided execution list from 'exec'.\n"); 
    if (have_outfile)
    fprintf(outfile,"Using the user provided execution list from 'exec'.\n"); 
    exec = parse_var(&nexec,MXEXEC,"EXEC");
    /* turn off auto-run of $input command */
    auto_input = 0;
  }
  else { /* Default sequence of modules from psi.dat */
    /* construct the name of the procedure from dertyp, reftyp, and wfn */
    strcpy(proced,"");
    if (direct == 1)
      strcat(proced,"DIRECT");
    if (strcmp(wfn,"CASSCF")==0 || strcmp(wfn,"RASSCF")==0)
      strcat(proced,"DETCAS");
    else
      strcat(proced,wfn);
    strcat(proced,reftyp);
    /*
      For some cases do not need to append DERTYPE since it's irrelevant
      and determined by JobType
    */
    if (JobType != DBOC) {
      if (strcmp(dertyp,"NONE")==0) strcat(proced,"ENERGY");
      else strcat(proced,dertyp);
    }
    
    switch (JobType) {
      case SP:
        jobtype = strdup("SP"); break;
      case OPT:
        jobtype = strdup("OPT"); break;
      case FREQ:
        jobtype = strdup("FREQ"); break;
      case SYMM_FC:
        jobtype = strdup("SYMM_FC"); break;
      case FC:
        jobtype = strdup("FC"); break;
      case DISP:
        jobtype = strdup("DISP"); break;
      case OEPROP:
        jobtype = strdup("OEPROP"); break;
      case DBOC:
        jobtype = strdup("DBOC"); break;
    }

    /* Unless jobtype = SP, append it to the calculation type */    
    if (strcmp(jobtype,"SP"))
      strcat(proced,jobtype);

    if (auto_check || check) {
      fprintf(stdout, "Calculation type string = %s\n", proced);
      if (have_outfile)
      fprintf(outfile, "Calculation type string = %s\n", proced);
      fprintf(stdout, "Wavefunction            = %s\n", wfn);
      if (have_outfile)
      fprintf(outfile, "Wavefunction            = %s\n", wfn);
      fprintf(stdout, "Reference               = %s\n", reftyp);
      if (have_outfile)
      fprintf(outfile, "Reference               = %s\n", reftyp);
      fprintf(stdout, "Jobtype                 = %s\n", jobtype);
      if (have_outfile)
      fprintf(outfile, "Jobtype                 = %s\n", jobtype);
      fprintf(stdout, "Dertype                 = %s\n", dertyp);
      if (have_outfile)
      fprintf(outfile, "Dertype                 = %s\n", dertyp);
      fprintf(stdout, "Direct                  = %s\n", 
        (direct ? "true" : "false"));
      if (have_outfile)
      fprintf(outfile, "Direct                  = %s\n", 
        (direct ? "true" : "false"));
    }

    if (!ip_exist(proced,0)) {
      fprintf(stdout,"Error: Did not find a valid calculation type,\n");
      fprintf(stdout,
              "probably because of incompatible wfn,dertype,reference.\n");
      fprintf(stdout,"Full calculation type is: %s\n", proced); 
      psi3_abort();
    }
    exec = parse_var(&nexec,MXEXEC,proced);
    /* fprintf(outfile, "nexec = %d\n", nexec); */
  } /* end default sequence from psi.dat */


  if (auto_check || check) {
    fprintf(stdout,"\n'CHECK' is YES, so nothing will be executed.\n");
    fprintf(stdout,
      "\nThe following programs would otherwise be executed:\n\n");
    if (have_outfile) {
    fprintf(outfile,"\n'CHECK' is YES, so nothing will be executed.\n");
    fprintf(outfile,
      "\nThe following programs would otherwise be executed:\n\n");
    }
  }
  else {
    fprintf(stdout,"\nThe following programs will be executed:\n\n");
    if (have_outfile)
      fprintf(outfile,"\nThe following programs will be executed:\n\n");
  }

  if(auto_input) {
    /* set up the "input" program execution, which should occur before the 
       rest of the procedure */
    input_exec = parse_var(&nexec_input, MXEXEC, "INPUT");
    fprintf(stdout, " %s\n", input_exec[0]);
    if (have_outfile)
    fprintf(outfile, " %s\n", input_exec[0]);
  }

  rdepth = 0;
  for (i=0; i<nexec; i++) {
    if (strcmp(exec[i],"END")==0) rdepth=rdepth-1;
    for (j=0; j<rdepth; j++) {
      fprintf(stdout," ");
      if (have_outfile) fprintf(outfile," ");
    }
    if (strcmp(exec[i],"REPEAT")==0) {
      rdepth=rdepth+1;
      fprintf(stdout," %s %s\n",exec[i],exec[i+1]);
      if (have_outfile) fprintf(outfile," %s %s\n",exec[i],exec[i+1]);
      i++;
    }
    else if (strlen(exec[i]) != 0)
      fprintf(stdout," %s\n",exec[i]);
      if (have_outfile) fprintf(outfile," %s\n",exec[i]);
  }

  fprintf(stdout,"\n");
  if (have_outfile) fprintf(outfile,"\n");
  
  /* don't write in between, b/c modules will be writing there...*/
  if (have_outfile) fclose(outfile);
  
  if(auto_input && !check && !auto_check) execut(input_exec,1,0);
  if (!check && !auto_check) execut(exec,nexec,0); 

  /* clean up and free memory */
  free(wfn);
  free(dertyp);
  free(reftyp);
  if (!ip_exist("EXEC",0)) free(jobtype);
  for (i=0; i<nexec; i++) free(exec[i]);
  ip_done();

  /* Normal completion */
  if (have_outfile) {
    outfile = fopen(ofname, "a");
    if (outfile == NULL) {
      fprintf(stderr, "Error: could not open output file %s\n", ofname);
      exit(1);
    }

    fprintf(outfile,"\n                          --------------------------\n");
    fprintf(outfile,"                          PSI3 Computation Completed\n");
    fprintf(outfile,"                          --------------------------\n");
    tstop_psi3driver(outfile);
  }

  exit(0);
}


namespace psi { namespace psi3 {

void psi3_abort(void)
{
  if (outfile)
    fprintf(outfile,"\nPSI3 exiting.\n");
  else
    fprintf(stderr, "\nPSI3 exiting.\n");

  abort();
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
      fprintf(stdout,"%s%s\n", spaces, exec[i]);
      runcmd(&errcod,exec[i]);

      /* fprintf(stderr,"%serrcod before filter is %d\n",spaces,errcod); */

      /* exited with signal */
      if (WIFSIGNALED(errcod)) {
        fprintf(stdout,"\nCommand %s was terminated with signal %d",
                exec[i], WTERMSIG(errcod));
	psi3_abort();
      }

      errcod = WEXITSTATUS(errcod);

      /* fprintf(stderr,"%serrcod after filter is %d\n",spaces,errcod); */

      /* if we're getting ndisp from optking */
      if ( !(strncmp("optking --disp_irrep",exec[i],
           strlen("optking --disp_irrep")))
        || !(strcmp(exec[i],"optking --disp_nosymm"))
        || !(strcmp(exec[i],"optking --disp_freq_grad_cart")) 
        || !(strcmp(exec[i],"optking --disp_freq_energy_cart")) 
        ) {
        /* fprintf(stdout, "optking got exit code %d\n", errcod); */
        for (j=i; j<nexec; j++) {
          /* fprintf(stdout,"Scanning exec[%d] = %s\n", j, exec[j]);  */
          if (strcmp(exec[j],"NUM_DISP")==0) {	
            sprintf(exec[j],"%d",errcod);
         /*   sprintf(exec[j], "%d", get_ndisp() ); */
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
          fprintf(stdout,"PSI3 caught ENDLOOP signal but no END below!\n");
          psi3_abort();
        }
        else return(-j);
      }
      /* if we got a 'command failed' failed flag */
      else if (errcod == PSI_RETURN_FAILURE) {
        fprintf(stdout,"\nCommand %s has returned a fail status.", exec[i]);
	psi3_abort();
      }
    }
    else {
      sscanf(exec[i+1], "%d", &nrep);
      fprintf(stdout, "%sREPEAT (%d)\n", spaces, nrep);
      for (j=0; j<nrep; j++) {
        fprintf(stdout, "%s CYCLE %d\n", spaces, j+1);
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

/* this code is essentially the same as in the 
psi_start(&infile,&outfile,&psi_file_prefix,) lib funct 
*/
int parse_cmdline(int argc, char *argv[])
{
  int i,errcod;
  int found_if_np = 0;          /* found input file name without -i */
  int found_of_np = 0;          /* found output file name without -o */
  int found_fp_np = 0;          /* found file prefix name without -p */
  int found_if_p = 0;           /* found input file name with -i */
  int found_of_p = 0;           /* found output file name with -o */
  int found_fp_p = 0;           /* found file prefix name with -p */
  const char *ifname=NULL, *fprefix=NULL; 
  char *arg;

  /* defaults */
  auto_input = 1;
  auto_check = 0;
  called_from_dboc = 0;
  call_done = 1;
  have_outfile = 1;

  /* process command-line arguments in sequence */
  for(i=1; i<argc; i++) {
    arg = argv[i];

    if(!strcmp(arg,"--noinput") || !strcmp(arg,"-n")) {
      auto_input = 0;
    }
    else if(!strcmp(arg,"--check") || !strcmp(arg,"-c")) {
      auto_check = 1;
    }
    /* check if called by dboc, which means ignore jobtype=dboc keyword
       and pretend jobtype=sp */
    else if(!strcmp(arg,"--dboc") || !strcmp(arg,"-d")) {
      called_from_dboc = 1;
    }
    /* check if $done command should not be executed */
    else if(!strcmp(arg,"--messy") || !strcmp(arg,"-m")) {
      call_done = 0;
    }
    else if ((!strcmp(arg,"-f") || !strcmp(arg,"-i")) && !found_if_p) {
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
    else if (!strcmp(arg,"-rp") || !strcmp(arg,"--randomprefix")) {

      char rand_prefix[8];
      pid_t pid;

      pid = getpid();

      errcod = sprintf(rand_prefix,"psi%05d",pid);

      setenv("PSI_PREFIX",rand_prefix,1);
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

  /* if some arguments were not specified on command-line - 
     check the environment */
  if (ifname == NULL)
    ifname = getenv("PSI_INPUT");
  if (ofname == NULL)
    ofname = getenv("PSI_OUTPUT");
  if (fprefix == NULL)
    fprefix = getenv("PSI_PREFIX");
 
  /* set the environmental variables the modules will look for */ 
  if (ifname != NULL) {
#if HAVE_PUTENV
    sprintf(tmpstr_input,"PSI_INPUT=%s",ifname);
    putenv(tmpstr_input);
#elif HAVE_SETENV
    setenv("PSI_INPUT",ifname,1);
#else
#error "Have neither putenv nor setenv. Something must be very broken on this system."
#endif
    infile = fopen(ifname,"r");
  }
  else {
    infile = fopen("input.dat","r");
    ifname = "input.dat";
  }
  if (infile == NULL) {
    fprintf(stderr, "Error: could not open input file %s\n",ifname);
    return(0);
  }
  if (ofname != NULL) {
#if HAVE_PUTENV
    sprintf(tmpstr_output,"PSI_OUTPUT=%s",ofname);
    putenv(tmpstr_output);
#elif HAVE_SETENV
    setenv("PSI_OUTPUT",ofname,1);
#else
#error "Have neither putenv nor setenv. Something must be very broken on this system."
#endif
  }
  else {
    ofname = "output.dat";
  }

  // CDS 3/08
  //outfile = stdout;
  if ((ofname[0]=='-' && ofname[1]=='\x0') || auto_check==1) {
    outfile=stdout;
    have_outfile = 0;
  }
  else {
    outfile = fopen(ofname, "w");
    if (outfile == NULL) {
      fprintf(stderr, "Error: could not open output file %s\n", ofname);
      return(0);
    }
  }

  if (fprefix != NULL) {
#if HAVE_PUTENV
    sprintf(tmpstr_prefix,"PSI_PREFIX=%s",fprefix);
    putenv(tmpstr_prefix);
#elif HAVE_SETENV
    setenv("PSI_PREFIX",fprefix,1);
#else
#error "Have neither putenv nor setenv. Something must be very broken on this system."
#endif
  }

  /*RAK */
  /* this hack works but doesn't pick up psi prefixes
  fclose(infile);
  psi_start(&infile,&outfile,&psi_file_prefix,0,NULL,0);
  ip_done();
  fclose(outfile);
  fclose(infile);
  if (ifname != NULL)
    infile = fopen(ifname,"r");
  else
    infile = fopen("input.dat","r");
  outfile = stdout;
  */
  return(1);
}

int get_ndisp(void) {
  int ndisp;
  FILE *outfile_psi3;

  outfile_psi3 = outfile;

  /* psi_file_prefix = getenv("PSI_PREFIX");
  printf("%s\n", psi_file_prefix);
  psi_start(&infile,&outfile,&psi_file_prefix,0,NULL,0); */

  /* need to remove psi_start */
  psio_init(); psio_ipv1_config();
  psio_open(PSIF_OPTKING, PSIO_OPEN_OLD);
  psio_read_entry(PSIF_OPTKING, "OPT: Num. of disp.", (char *) &(ndisp), 
    sizeof(int));
  psio_close(PSIF_OPTKING,1);
  psio_done();
  printf("ndisp: %d\n",ndisp);
  outfile = outfile_psi3;
  return ndisp;
}

void print_welcome(FILE *outfile)
{
fprintf(outfile, 
"    -----------------------------------------------------------------------    \n");
fprintf(outfile,
"            PSI3: An Open-Source Ab Initio Electronic Structure Package \n");
fprintf(outfile,
"                                Version 3.4 Alpha\n\n");
fprintf(outfile,
"    T. D. Crawford, C. D. Sherrill, E. F. Valeev, J. T. Fermann, R. A. King,\n");
fprintf(outfile,
"    M. L. Leininger, S. T. Brown, C. L. Janssen, E. T. Seidl, J. P. Kenny,\n");
fprintf(outfile,
"    and W. D. Allen, J. Comput. Chem. 28, 1610-1616 (2007)\n");
fprintf(outfile,
"    -----------------------------------------------------------------------    \n");
}

void tstart_psi3driver(FILE *outfile)
{
  using namespace psi::psi3;

  int i,error;
  char *name;
  name = (char *) malloc(40 * sizeof(char));
  error = gethostname(name, 40);
  if(error != 0) strncpy(name,"nohostname", 11);

  time_start = time(NULL);

  for (i=0; i < 78 ; i++) {
    fprintf(outfile,"*");
  }
  fprintf(outfile,"\n");

  fprintf(outfile,"PSI3 started on %s at %s\n", name, ctime(&time_start));

  free(name);
}

void tstop_psi3driver(FILE *outfile)
{
  using namespace psi::psi3;

  int i;
  int error;
  time_t total_time;
  struct tms total_tmstime;
  char *name;
  double user_s, sys_s;

  name = (char *) malloc(40 * sizeof(char));
  error = gethostname(name, 40);
  if(error != 0) strncpy(name,"nohostname", 11);

  time_end = time(NULL);
  total_time = time_end - time_start;

  times(&total_tmstime);
  const long clk_tck = sysconf(_SC_CLK_TCK);
  user_s = ((double) total_tmstime.tms_utime)/clk_tck;
  sys_s = ((double) total_tmstime.tms_stime)/clk_tck;

  fprintf(outfile,"\n");
  fprintf(outfile,"PSI3 stopped on %s at %s\n", name, ctime(&time_end));
  fprintf(outfile,"Total PSI3 wall time %10d seconds = %10.2f minutes\n",
          total_time, ((double) total_time)/60.0);

  for (i=0; i < 78 ; i++) {
    fprintf(outfile,"*");
  }

  free(name);
}

}} // namespace psi::psi3
