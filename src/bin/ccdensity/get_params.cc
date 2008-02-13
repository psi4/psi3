/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void get_params()
{
  int errcod, tol;
  char *junk;

  errcod = ip_string("WFN", &(params.wfn), 0);

  if(!strcmp(params.wfn,"CC2") || !strcmp(params.wfn,"EOM_CC2")) {
    psio_read_entry(CC_INFO, "CC2 Energy", (char *) &(moinfo.ecc),
                    sizeof(double));
    fprintf(outfile,  "\tCC2 energy          (CC_INFO) = %20.15f\n",moinfo.ecc);
    fprintf(outfile,  "\tTotal CC2 energy    (CC_INFO) = %20.15f\n",
            moinfo.eref+moinfo.ecc);
  }
  else if(!strcmp(params.wfn,"CCSD") || !strcmp(params.wfn,"EOM_CCSD")) {
    psio_read_entry(CC_INFO, "CCSD Energy", (char *) &(moinfo.ecc),
                    sizeof(double));
    fprintf(outfile,  "\tCCSD energy         (CC_INFO) = %20.15f\n",moinfo.ecc);
    fprintf(outfile,  "\tTotal CCSD energy   (CC_INFO) = %20.15f\n",
            moinfo.eref+moinfo.ecc);
  }
  else if(!strcmp(params.wfn,"CC3") || !strcmp(params.wfn,"EOM_CC3")) {
    psio_read_entry(CC_INFO, "CC3 Energy", (char *) &(moinfo.ecc),
                    sizeof(double));
    fprintf(outfile,  "\tCC3 energy          (CC_INFO) = %20.15f\n",moinfo.ecc);
    fprintf(outfile,  "\tTotal CC3 energy    (CC_INFO) = %20.15f\n",
            moinfo.eref+moinfo.ecc);
  }

  fflush(outfile);

  params.tolerance = 1e-14;
  errcod = ip_data("TOLERANCE","%d",&(tol),0);
  if(errcod == IPE_OK) params.tolerance = 1.0*pow(10.0,(double) -tol);

  fndcor(&(params.memory),infile,outfile);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);

  params.aobasis = 0;
  errcod = ip_boolean("AO_BASIS", &(params.aobasis),0);

  params.ael = 0;
  errcod = ip_boolean("AEL", &(params.ael),0);

  if(ip_exist("GAUGE",0)) {
    ip_string("GAUGE",&(params.gauge), 0);
    if(strcmp(params.gauge,"LENGTH") && strcmp(params.gauge,"VELOCITY")) {
      printf("Invalid choice of gauge: %s\n", params.gauge);
      exit(PSI_RETURN_FAILURE);
    }
  }
  else params.gauge = strdup("LENGTH");
  
  /*** determine DERTYPE from input */
  params.dertype = 0;
  if(ip_exist("DERTYPE",0)) {
    errcod = ip_string("DERTYPE", &(junk),0);
    if(!strcmp(junk,"NONE")) params.dertype = 0;
    else if(!strcmp(junk,"FIRST")) params.dertype = 1;
    else if(!strcmp(junk,"RESPONSE")) params.dertype = 3; 
    else {
      printf("Invalid value of input keyword DERTYPE: %s\n", junk);
      exit(PSI_RETURN_FAILURE);
    }
    free(junk);
  }
  else { /* infer DERTYPE from JOBTYPE */
    errcod = ip_string("JOBTYPE", &(junk),0);
    if ( !strcmp(junk,"SP") ) params.dertype = 0;
    else if ( !strcmp(junk, "OEPROP") ) params.dertype = 0;
    else if ( !strcmp(junk, "OPT") ) params.dertype = 1;
    else {
      printf("Not sure what to do with missing DERTYPE and this JOBTYPE: %s\n", junk);
      exit(PSI_RETURN_FAILURE);
    }
    free(junk);
  }

  if ( (params.dertype == 1) || (!strcmp(params.wfn, "CCSD_MVD")) )
    params.relax_opdm = 1;  /* default for gradients, or MVD correction */
  else
    params.relax_opdm = 0;  /* otherwise, default is relax_opdm off */

  if(params.transition) 
    params.relax_opdm = 0;

  errcod = ip_boolean("RELAX_OPDM", &(params.relax_opdm),0);
  if ( (params.onepdm) && (params.relax_opdm) ) { /* can't do relaxation without twopdm */
    fprintf(outfile,"\tTurning orbital relaxation off since only onepdm is requested.\n");
    params.relax_opdm = 0;
  }

  if ( (!strcmp(params.wfn,"EOM_CCSD")) && (params.dertype == 0) )
    params.connect_xi = 0;
  else
    params.connect_xi = 1;
  errcod = ip_boolean("CONNECT_XI",&(params.connect_xi),0);

  
  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tTolerance        = %3.1e\n", params.tolerance);
  fprintf(outfile, "\tCache Level      = %1d\n", params.cachelev);
  fprintf(outfile, "\tAO Basis         = %s\n", 
          params.aobasis ? "Yes" : "No");
  fprintf(outfile, "\tOPDM Only        = %s\n", 
	  params.onepdm ? "Yes" : "No");
  fprintf(outfile, "\tRelax OPDM       = %s\n", 
          params.relax_opdm ? "Yes" : "No");
  fprintf(outfile, "\tCompute Xi       = %s\n", 
          (params.calc_xi) ? "Yes" : "No");
  fprintf(outfile, "\tUse Zeta         = %s\n", 
          (params.use_zeta) ? "Yes" : "No");
  fprintf(outfile, "\tXi connected     = %s\n", 
          (params.connect_xi) ? "Yes" : "No");
  fflush(outfile);
}


}} // namespace psi::ccdensity
