/*******************************************************************
 *  PSIMRCC (2008) by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ******************************************************************/

/**
 *  @defgroup PSIMRCC PSIMRCC is a code for SR/MRCC computations
 *  @file psimrcc.cpp
 *  @ingroup (PSIMRCC)
 *  @brief Contains main() and global variables
*/

// Standard libraries 
#include <iostream>
#include <complex>
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
#include <libqt/qt.h>

#include "blas.h"
#include "debugging.h"
#include "git.h"
#include "main.h"
#include "sort.h"
#include "mrcc.h"
#include "memory_manager.h"
#include "psimrcc.h"
#include "transform.h"

// PSI FILES
FILE  *infile, *outfile;
char  *psi_file_prefix;
extern "C" {
  char* gprgid();
}

using namespace std;

MOInfo              *moinfo;

namespace psi{ namespace psimrcc{

// Global variables
Timer               *global_timer;
Debugging           *debugging;
MemoryManager       *mem;
CCBLAS              *blas;
CCSort              *sorter;
CCTransform         *trans = NULL;

void read_calculation_options();
}} /* End Namespaces */

/**
 * The main function
 * @param argc
 * @param argv[]
 * @return PSI_RETURN_SUCCESS if the program ran without any problem
 */
int main(int argc, char *argv[])
{
  using namespace psi::psimrcc;
  init_psi(argc,argv);

  global_timer = new Timer;
  read_calculation_options();

  debugging = new Debugging;

  mem    = new MemoryManager();
  moinfo = new MOInfo();
  moinfo->setup_model_space();

  run_psimrcc();

  fprintf(outfile,"\n\n\tPSIMRCC Execution Ended.");
  fprintf(outfile,"\n\tWall Time = %20.6f s",global_timer->get());
  fprintf(outfile,"\n\tGEMM Time = %20.6f s",moinfo->get_dgemm_timing());
  fflush(outfile);

  mem->MemCheck(outfile);

  delete moinfo;
  delete mem;
  delete global_timer;
  close_psi();
  return PSI_RETURN_SUCCESS;
}


namespace psi{ namespace psimrcc{

void read_calculation_options()
{
  options_add_int("CORR_CHARGE",0);
  options_add_int("DEBUG",0);
  options_add_int("DAMPING_FACTOR",0);
  options_add_int("MAXDIIS",7);
  options_add_int("MEMORY",1800);
  options_add_int("NUM_THREADS",1);
  options_add_int("NEL",0);
  options_add_int("ROOT",1);
  options_add_int("E_CONVERGENCE",9);
  options_add_int("PT_E_CONVERGENCE",9);
  options_add_int("MAX_ITERATIONS",100);
  options_add_int("DENOMINATOR_SHIFT",0);
  options_add_int("START_DIIS",2);

  options_add_bool("DIIS_TRIPLES",false);
  options_add_bool("LOCK_SINGLET",false);
  options_add_bool("MP2_GUESS",true);
  options_add_bool("ONLY_CLOSED_SHELL",false);
  options_add_bool("USE_DIIS",true);
  options_add_bool("USE_SPIN_SYMMETRY",true);
  options_add_bool("ZERO_INTERNAL_AMPS",true);
  options_add_bool("COUPLING_TERMS",true);
  options_add_bool("PRINT_HEFF",false);

  options_add_str_with_choices("PT_ENERGY","SECOND_ORDER","SECOND_ORDER SCS_SECOND_ORDER PSEUDO_SECOND_ORDER SCS_PSEUDO_SECOND_ORDER");
  options_add_str_with_choices("CORR_WFN","CCSD","PT2 CCSD CCSD_T CCSDT-1A CCSDT-1B CCSDT-2 CCSDT-3 CCSDT MP2-CCSD");
  options_add_str_with_choices("CORR_REFERENCE","GENERAL","RHF ROHF TCSCF MCSCF GENERAL");
  options_add_str_with_choices("CORR_ANSATZ","MK","SR MK BW APBW");
  options_add_str_with_choices("COUPLING","CUBIC","NONE LINEAR QUADRATIC CUBIC");
  options_add_str_with_choices("WFN_SYM","1","A AG AU AP APP A1 A2 B BG BU B1 B2 B3 B1G B2G B3G B1U B2U B3U 0 1 2 3 4 5 6 7 8");

  // MP2_CCSD
  options_add_str_with_choices("MP2_CCSD_METHOD","II","I IA II");

  options_read();
  options_print();
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
  ip_cwk_add(const_cast<char*>(":MRCC"));
  ip_cwk_add(const_cast<char*>(":PSIMRCC"));
  ip_cwk_add(const_cast<char*>(":PSI"));
  ip_cwk_add(const_cast<char*>(":INPUT"));

  
  options_init();

  tstart(outfile);
  
  psio_open(PSIF_PSIMRCC_INTEGRALS,PSIO_OPEN_NEW);

  

  fprintf(outfile,"\n  MRCC          MRCC");
  fprintf(outfile,"\n   MRCC  MRCC  MRCC");
  fprintf(outfile,"\n   MRCC  MRCC  MRCC      PSIMRCC Version 0.8.0, August, 2008");
  fprintf(outfile,"\n   MRCC  MRCC  MRCC      Multireference Coupled Cluster, written by");
  fprintf(outfile,"\n     MRCCMRCCMRCC        Francesco A. Evangelista and Andrew C. Simmonett");
  fprintf(outfile,"\n         MRCC            Compiled on %s at %s",__DATE__,__TIME__);
  fprintf(outfile,"\n         MRCC            id =%s",GIT_ID);
  fprintf(outfile,"\n       MRCCMRCC");          
}

/**
 * Close psi by calling psio_done() and psi_stop()
 */
void close_psi()
{
  /***********************
    Close the checkpoint
  ***********************/
  fflush(outfile);
  options_close();
  chkpt_close();
  psio_close(PSIF_PSIMRCC_INTEGRALS,1);
  psio_done();
  tstop(outfile);
  psi_stop(infile,outfile,psi_file_prefix);
}

}} /* End Namespaces */

/**
 * @return program ID
 */
char* gprgid()
{
  return(const_cast<char*>("PSIMRCC"));
}
