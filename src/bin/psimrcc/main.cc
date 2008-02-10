/***************************************************************************
 *  PSIMRCC
 *  Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

/**
 *  @defgroup PSIMRCC PSIMRCC is a code for SR/MRCC computations
 *  @file psimrcc.cpp
 *  @ingroup (PSIMRCC)
 *  @brief Contains main() and global variables
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <complex>


#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>

#include "git.h"
#include "main.h"
#include "psimrcc.h"
#include "sq.h"
#include "blas.h"
#include "sort.h"
#include "mrcc.h"
#include "memory_manager.h"
#include "transform.h"
#include "moinfo.h"
#include "calculation_options.h"
#include "debugging.h"

#include "utilities.h"


// PSI FILES
FILE  *infile, *outfile;
char  *psi_file_prefix;
extern "C" {
  char* gprgid();
}

using namespace std;

namespace psi{ namespace psimrcc{

// Global variables
Timer               *global_timer;
CalculationOptions  *options;
Debugging           *debugging;
MemoryManager       *mem;
MOInfo              *moinfo;
CCBLAS              *blas;
CCSort              *sorter;
CCTransform         *trans = NULL;

void read_calculation_options();
}} /* End Namespaces */

/**
 * The main function
 * @param argc
 * @param argv[]
 * @return EXIT_SUCCESS if the program ran without any problem
 */
int main(int argc, char *argv[])
{
  using namespace psi::psimrcc;
  init_psi(argc,argv);

  global_timer = new Timer;
  options      = new CalculationOptions;
  read_calculation_options();

  debugging = new Debugging;

  mem    = new MemoryManager();
  moinfo = new MOInfo(argc,argv,":MRCC");
  moinfo->setup_model_space();

  if(options->get_str_option("CORR_CODE")=="PSIMRCC")
    run_psimrcc();
  if(options->get_str_option("CORR_CODE")=="SQ")
    run_sq();
//   if(options->get_str_option("CORR_CODE")=="MCSCF")
//     run_tcscf();
//   switch (options->get_str_option("CORR_CODE")){
//     case "PSIMRCC" :
//       run_psimrcc();
//       break;
//     case "SQ" : 
//       run_sq();
//       break;
//     case "MCSCF" : 
//       break;
//     default:
//       break;
//   }

  fprintf(outfile,"\n\n\tPSIMRCC Execution Ended.");
  fprintf(outfile,"\n\tWall Time = %20.6f s",global_timer->get());
  fprintf(outfile,"\n\tGEMM Time = %20.6f s",moinfo->get_dgemm_timing());
  fflush(outfile);

  double* my_vector;
  allocate1(double,my_vector,10);

  release1(my_vector);

  mem->MemCheck(outfile);



  delete moinfo;
  delete mem;
  delete options;
  delete global_timer;
  close_psi();
  return EXIT_SUCCESS;
}


namespace psi{ namespace psimrcc{

void read_calculation_options()
{
  options->add_int_option("CORR_CHARGE",0);
  options->add_int_option("DEBUG",0);
  options->add_int_option("DAMPING_FACTOR",0);
  options->add_int_option("MAXDIIS",7);
  options->add_int_option("MEMORY",1800);
  options->add_int_option("MULTP",1);
  options->add_int_option("NUM_THREADS",1);
  options->add_int_option("NEL",0);
  options->add_int_option("ROOT",1);
  options->add_int_option("E_CONVERGENCE",9);
  options->add_int_option("PT_E_CONVERGENCE",9);
  options->add_int_option("MAX_ITERATIONS",100);
  options->add_int_option("DENOMINATOR_SHIFT",0);

  options->add_bool_option("DIIS_TRIPLES",false);
  options->add_bool_option("LOCK_SINGLET",false);
  options->add_bool_option("MP2_GUESS",true);
  options->add_bool_option("ONLY_CLOSED_SHELL",false);
  options->add_bool_option("USE_DIIS",true);
  options->add_bool_option("USE_SPIN_SYMMETRY",true);
  options->add_bool_option("ZERO_INTERNAL_AMPS",true);

  options->add_str_option_with_choices("PT_ENERGY","SECOND_ORDER","SECOND_ORDER SCS_SECOND_ORDER PSEUDO_SECOND_ORDER SCS_PSEUDO_SECOND_ORDER");
  options->add_str_option_with_choices("CORR_WFN","CCSD","PT2 CCSD CCSD_T CCSDT-1A CCSDT-1B CCSDT-2 CCSDT-3 CCSDT");
  options->add_str_option_with_choices("CORR_ANSATZ","MK","SR MK BW APBW");
  options->add_str_option_with_choices("CORR_CODE","PSIMRCC","PSIMRCC SQ LOOPMRCC");
  options->add_str_option_with_choices("COUPLING","CUBIC","NONE LINEAR QUADRATIC CUBIC");
  options->add_str_option_with_choices("WFN_SYM","1","A AG AU AP APP A1 A2 B BG BU B1 B2 B3 B1G B2G B3G B1U B2U B3U 0 1 2 3 4 5 6 7 8");
  options->read_options();
  options->print();
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
  ip_cwk_add(":PSI");
  ip_cwk_add(":INPUT");
  ip_cwk_add(":TRANSQT");
  ip_cwk_add(":MRCC");
  psio_open(MRCC_ON_DISK,PSIO_OPEN_NEW);


  fprintf(outfile,"\n  MRCC          MRCC");
  fprintf(outfile,"\n   MRCC  MRCC  MRCC");
  fprintf(outfile,"\n   MRCC  MRCC  MRCC      PSIMRCC Version 0.7.2, November 07, 2007");
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
  chkpt_close();
  psio_close(MRCC_ON_DISK,1);
  psio_done();
  psi_stop(infile,outfile,psi_file_prefix);
}

}} /* End Namespaces */

/**
 * @return program ID
 */
char* gprgid()
{
  return(":MRCC");
}
