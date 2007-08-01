
#include <stdlib.h>
#include <string.h>
#include <psifiles.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>

extern char* gprgid();

/*!
  This function initializes the libpsio library using the .psirc file + input file
*/
int psio_ipv1_config()
{
  int i, unit, nkwds;
  char* ip_token;
  char *name;
  char* userhome;
  FILE* psirc;
  char* filename;
  int errcod;

  _psio_error_exit_code_ = PSI_RETURN_FAILURE;

  /* Open user's general .psirc file, if exists */
  userhome = getenv("HOME");
  filename = (char*) malloc((strlen(userhome)+8)*sizeof(char));
  sprintf(filename, "%s%s", userhome, "/.psirc");
  psirc = fopen(filename, "r");
  if(psirc != NULL) {
    ip_append(psirc, stdout);
    fclose(psirc);
  }
  free(filename);

  /*
    transfer the necessary keywords from IPV1 to PSIO
  */
  /* need "NAME", "NVOLUME", and "VOLUMEX" -- total of 2+PSIO_MAXVOL keywords */
  nkwds = 2+PSIO_MAXVOL;
  char** kwds = (char**) malloc(nkwds*sizeof(char*));
  kwds[0] = strdup("NAME");
  kwds[1] = strdup("NVOLUME");
  for(i=1; i<=PSIO_MAXVOL; ++i) {
    char kwd[20];
    sprintf(kwd,"VOLUME%u",i);
    kwds[2+i] = strdup(kwd);
  }

  /* allocate ip_token
     conservative estimate for its length = strlen(gprgid())+80
  */
  ip_token = (char*) malloc( (strlen(gprgid())+80)*sizeof(char) );

  for(i=0; i<nkwds; ++i) {
    const char* kwd = kwds[i];

    /* unit and program specific */
    for(unit=0; unit<PSIO_MAXUNIT; ++unit) {
      sprintf(ip_token,":%s:FILES:FILE%u:%s",gprgid(),unit,kwd);
      errcod = ip_data(ip_token,"%s",name,0);
      if (errcod == IPE_OK) {
	psio_set_filescfg_kwd(gprgid(),kwd,unit,name);
	free(name);
      }
    }
    
    /* program specific */
    sprintf(ip_token,":%s:FILES:DEFAULT:%s",gprgid(),kwd);
    errcod = ip_data(ip_token,"%s",name,0);
    if (errcod == IPE_OK) {
      psio_set_filescfg_kwd(gprgid(),kwd,-1,name);
      free(name);
    }

    /* unit specific in PSI section */
    for(unit=0; unit<PSIO_MAXUNIT; ++unit) {
      sprintf(ip_token,":PSI:FILES:FILE%u:%s",unit,kwd);
      errcod = ip_data(ip_token,"%s",name,0);
      if (errcod == IPE_OK) {
	psio_set_filescfg_kwd("PSI",kwd,unit,name);
	free(name);
      }
    }

    /* in PSI section */
    sprintf(ip_token,":PSI:FILES:DEFAULT:%s",kwd);
    errcod = ip_data(ip_token,"%s",name,0);
    if (errcod == IPE_OK) {
      psio_set_filescfg_kwd("PSI",kwd,-1,name);
      free(name);
    }

    /* unit specific in DEFAULT section */
    for(unit=0; unit<PSIO_MAXUNIT; ++unit) {
      sprintf(ip_token,":DEFAULT:FILES:FILE%u:%s",unit,kwd);
      errcod = ip_data(ip_token,"%s",name,0);
      if (errcod == IPE_OK) {
	psio_set_filescfg_kwd("DEFAULT",kwd,unit,name);
	free(name);
      }
    }

    /* in DEFAULT section */
    sprintf(ip_token,":DEFAULT:FILES:DEFAULT:%s",kwd);
    errcod = ip_data(ip_token,"%s",name,0);
    if (errcod == IPE_OK) {
      psio_set_filescfg_kwd("DEFAULT",kwd,-1,name);
      free(name);
    }

  }

  /*
    implement some default PSI3 behavior:
    1) checkpoint file should by default be in "./"
   */
  {
    for(i=1; i<=PSIO_MAXVOL; ++i) {
      char kwd[20];
      sprintf(kwd,"VOLUME%u",i);
      psio_set_filescfg_kwd("DEFAULT",kwd,PSIF_CHKPT,"./");
    }
  }


  for(i=0; i<nkwds; ++i) {
    free(kwds[i]);
  }
  free(kwds);
  free(ip_token);

  return 0;
}
