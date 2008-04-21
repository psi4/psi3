//  MCSCF (2008) by Francesco Evangelista - frank@ccc.uga.edu
//  Center for Computational Chemistry, University of Georgia

/**
 *  @defgroup MCSCF MCSCF is a code for SCF/MCSCF computations
 *  @ingroup MCSCF
 *  @brief Contains main() and global variables
*/

// Standard libraries 
#include <iostream>
#include <cstdlib>

// PSI libraries 
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <libchkpt/chkpt.h>
#include <libmoinfo/libmoinfo.h>
#include <liboptions/liboptions.h>
#include <libutil/libutil.h>


#include "mcscf.h"
#include "memory_manager.h"
#include "git.h"
#include "scf.h"

// PSI FILES
FILE  *infile, *outfile;
char  *psi_file_prefix;
extern "C" {
  char* gprgid();
}

using namespace std;

MOInfoSCF              *moinfo_scf;
MemoryManager          *mem;

namespace psi{ namespace mcscf{

// Global variables
// Timer               *global_timer;
// Debugging           *debugging;


void add_calculation_options();
}} /* End Namespaces */

/**
 * The main function
 * @param argc
 * @param argv[]
 * @return PSI_RETURN_SUCCESS if the program ran without any problem
 */
int main(int argc, char *argv[])
{
  using namespace psi::mcscf;
  init_psi(argc,argv);

  mem    = new MemoryManager();
  moinfo_scf = new MOInfoSCF();

  if(options_get_str("REFERENCE") == "RHF"  ||
     options_get_str("REFERENCE") == "ROHF" ||
     options_get_str("REFERENCE") == "UHF"  ||
     options_get_str("REFERENCE") == "TWOCON"){
    SCF scf;
    scf.compute_energy();
  }else if(options_get_str("REFERENCE") == "MCSCF"){
    fprintf(outfile,"\n\nREFERENCE = MCSCF not implemented yet");
    fflush(outfile);
    return PSI_RETURN_FAILURE;
  }

  if(options_get_int("DEBUG") > 0)
    mem->MemCheck(outfile);
  delete moinfo_scf;
  delete mem;
  close_psi();
  return PSI_RETURN_SUCCESS;
}


namespace psi{ namespace mcscf{

void add_calculation_options()
{

  options_add_int("CONVERGENCE",9);
  options_add_int("DAMPING_FACTOR",0);
  options_add_int("DEBUG",0);
  options_add_int("DENOMINATOR_SHIFT",0);
  options_add_int("MAX_ITERATIONS",100);
  options_add_int("MEMORY",1800);
  options_add_int("NDIIS",7);
  options_add_int("ROOT",1);
  options_add_int("START_FAVG",5);
  options_add_int("TURN_ON_ACTV",5);


  options_add_bool("CI_DIIS",true);
  options_add_bool("USE_DIIS",true);
  options_add_bool("READ_MOS",true);
  options_add_bool("USE_FAVG",false);


  options_add_str_with_choices("REFERENCE","RHF","RHF ROHF UHF TWOCON MCSCF GENERAL");
  options_add_str_with_choices("WFN_SYM","1","A AG AU AP APP A1 A2 B BG BU B1 B2 B3 B1G B2G B3G B1U B2U B3U 0 1 2 3 4 5 6 7 8");
}

/**
 * Start Psi3, draw the program logo, version, and compilation details
 * @param argc
 * @param argv[]
 */
void init_psi(int argc, char *argv[])
{
  int num_extra_args=0;
  char**  extra_args;

  extra_args = new char*[argc];

  for(int i=1; i<argc; i++) {
    extra_args[num_extra_args++] = argv[i];
/*  Template for argument parsing
    if(strcmp(argv[i], "--opdm") == 0) {
    }
    else {
      extra_args[num_extra_args++] = argv[i];
    }
*/
  }
  psi_start(&infile,&outfile,&psi_file_prefix,num_extra_args,extra_args,0);
  delete[] extra_args;

  psio_init();
  psio_ipv1_config();
  chkpt_init(PSIO_OPEN_OLD);
  ip_cwk_clear();
  ip_cwk_add(const_cast<char*>(":PSI"));
  ip_cwk_add(const_cast<char*>(":SCF"));
  ip_cwk_add(const_cast<char*>(":MCSCF"));

  fprintf(outfile,"\n  MCSCF Version 0.1.0, April, 2008");
  fprintf(outfile,"\n  Francesco Evangelista");
  fprintf(outfile,"\n  Compiled on %s at %s",__DATE__,__TIME__);
  fprintf(outfile,"\n  id =%s",GIT_ID); 


  options_init();
  add_calculation_options();
  options_read();
  options_print();

  psio_open(PSIF_MCSCF,PSIO_OPEN_NEW);
}

/**
 * Close psi by calling psio_done() and psi_stop()
 */
void close_psi()
{
  fprintf(outfile,"\n\n  MCSCF Execution Completed.");
  fflush(outfile);

  options_close();

  chkpt_close();

  psio_close(PSIF_MCSCF,1);

  psio_done();

  psi_stop(infile,outfile,psi_file_prefix);
}

}} /* End Namespaces */

/**
 * @return program ID
 */
char* gprgid()
{
  return(const_cast<char*>(":MCSCF"));
}
