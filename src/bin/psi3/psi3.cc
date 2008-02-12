/*! \file psi3.cc
    \ingroup PSI3
    \brief Enter brief description of file here 
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <psifiles.h>          /* where return values are */
#include <sys/types.h>
#include <sys/wait.h>
#include <libpsio/psio.h>
#include <psifiles.h>

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
int get_ndisp(void);
/* 
  the following must stay in scope throughout the run, else the 
  environmental variables will vanish!
*/

char tmpstr_input[200], tmpstr_output[200], tmpstr_prefix[200];

void psi3_abort(void);
int execut(char **module_names, int num_modules, int depth);
extern char **parse_var(int *nvars, int mxvars, char *name);
extern void runcmd(int *errcod, char *cmd);
int parse_cmdline(int argc, char *argv[]);

/* boolean for automatically running input program first in all procedures */
int auto_input;
/* boolean for checking input w/o running any programs */
int auto_check;
/* boolean for whether called from dboc module (if yes - ignore jobtype=dboc) */
int called_from_dboc;
/* boolean for whether $done should be called */
int call_done;
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

  fprintf(outfile, "\n\n The PSI3 Execution Driver \n");

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
    fprintf(outfile, "Warning: didn't find psi.dat at %s\n", psidat_filename);
  }
  free(psidat_filename);

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
      fprintf(outfile, 
              "Error: A problem arose reading the optional keyword 'check'.\n");
      psi3_abort();
    }
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
    fprintf(outfile,"Error: Unrecognized calculation type %s\n", calctyp);
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
      fprintf(outfile,"Error: bad 'reference'.\n");
      fprintf(outfile,"Must be one of: RHF, ROHF, UHF, TWOCON\n");
      psi3_abort();
    }

  if ((strcmp(dertyp,"NONE")!=0) && (strcmp(dertyp,"FIRST")!=0) &&
      (strcmp(dertyp,"SECOND")!=0) && (strcmp(dertyp,"RESPONSE")!=0))
    {
      fprintf(outfile,"Error: bad 'dertype'.\n");
      fprintf(outfile,"Must be one of: NONE, FIRST, SECOND, or RESPONSE\n");
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
  else if (JobType == SYMM_FC) {
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
    else if (strcmp(dertyp,"RESPONSE")==0) fprintf(outfile,"response property");
    else fprintf(outfile, "unrecognized-dertype");
    fprintf(outfile, " computation.\n");
  }
  else if (JobType == OEPROP) {
    fprintf(outfile, "one-electron properties computation.\n");
  }
  else if (JobType == DBOC) {
    fprintf(outfile, 
      "Diagonal Born-Oppenheimer Correction (DBOC) computation.\n");
  }
  else { 
    fprintf(outfile, "calculation of unrecognized type.\n");
  }

  if (ip_exist("EXEC",0)) { /* User-specified EXEC statement */
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

    if (check || auto_check) {
      fprintf(outfile, "Calculation type string = %s\n", proced);
      fprintf(outfile, "Wavefunction            = %s\n", wfn);
      fprintf(outfile, "Reference               = %s\n", reftyp);
      fprintf(outfile, "Jobtype                 = %s\n", jobtype);
      fprintf(outfile, "Dertype                 = %s\n", dertyp);
      fprintf(outfile, "Direct                  = %s\n", 
        (direct ? "true" : "false"));
    }

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


  if (check || auto_check) {
    fprintf(outfile,"\n'CHECK' is YES, so nothing will be executed.\n");
    fprintf(outfile,"\nThe following programs would otherwise be executed:\n\n");
  }
  else
    fprintf(outfile,"\nThe following programs will be executed:\n\n");

  if(auto_input) {
    /* set up the "input" program execution, which should occur before the 
       rest of the procedure */
    input_exec = parse_var(&nexec_input, MXEXEC, "INPUT");
    fprintf(outfile, " %s\n", input_exec[0]);
  }

  rdepth = 0;
  for (i=0; i<nexec; i++) {
    if (strcmp(exec[i],"END")==0) rdepth=rdepth-1;
    for (j=0; j<rdepth; j++) fprintf(outfile," ");
    if (strcmp(exec[i],"REPEAT")==0) {
      rdepth=rdepth+1;
      fprintf(outfile," %s %s\n",exec[i],exec[i+1]);
      i++;
    }
    else if (strlen(exec[i]) != 0)
	fprintf(outfile," %s\n",exec[i]);
  }

  fprintf(outfile,"\n");
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
      fprintf(outfile,"%s%s\n", spaces, exec[i]);
      runcmd(&errcod,exec[i]);

      /* fprintf(stderr,"%serrcod before filter is %d\n",spaces,errcod); */

      /* exited with signal */
      if (WIFSIGNALED(errcod)) {
        fprintf(outfile,"\nCommand %s was terminated with signal %d",
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
        /* fprintf(outfile, "optking got exit code %d\n", errcod); */
        for (j=i; j<nexec; j++) {
          /* fprintf(outfile,"Scanning exec[%d] = %s\n", j, exec[j]);  */
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
          fprintf(outfile,"PSI3 caught ENDLOOP signal but no END below!\n");
          psi3_abort();
        }
        else return(-j);
      }
      /* if we got a 'command failed' failed flag */
      else if (errcod == PSI_RETURN_FAILURE) {
        fprintf(outfile,"\nCommand %s has returned a fail status.", exec[i]);
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

/* this code is essentially the same as in the psi_start(&infile,&outfile,&psi_file_prefix,) lib funct */
int parse_cmdline(int argc, char *argv[])
{
  int i,errcod;
  int found_if_np = 0;          /* found input file name without -i */
  int found_of_np = 0;          /* found output file name without -o */
  int found_fp_np = 0;          /* found file prefix name without -p */
  int found_if_p = 0;           /* found input file name with -i */
  int found_of_p = 0;           /* found output file name with -o */
  int found_fp_p = 0;           /* found file prefix name with -p */
  char *ifname=NULL, *ofname=NULL, *fprefix=NULL, *arg;

  /* defaults */
  auto_input = 1;
  auto_check = 0;
  called_from_dboc = 0;
  call_done = 1;

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
  outfile = stdout;
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

}} // namespace psi::psi3
