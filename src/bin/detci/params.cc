/*
** PARAMS.CC: File contains functions which get or print the running
**    parameters for the CI calculation.
**
** David Sherrill, 16 November 1994
**
*/

#define EXTERN

extern "C" {
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef AIX
#include <string.h>
#endif
#include <ip_libv1.h>
#include <libipv1/ip_data.gbl>
#include <libciomr.h>
#include <qt.h>
#include <file30.h>
#include <psifiles.h>
#include "structs.h"
#include "globals.h"
}


/*
** get_parameters(): Function gets the program running parameters such
**   as convergence.  These are stored in the Parameters data structure.
*/
void get_parameters(void)
{
   int i, errcod;
   int iopen=0, tval;
   char line1[133];
   
   ip_set_uppercase(1);
   ip_initialize(infile, outfile);
   ip_cwk_clear();
   ip_cwk_add(":DEFAULT");
   ip_cwk_add(":DETCI");


   /* default value of Ms0 depends on iopen but is modified below 
    * depending on value of opentype
    */
   file30_init();
   Parameters.Ms0 = !(file30_rd_iopen());
   file30_close();

   /* need to figure out wheter to filter tei's */
   errcod = ip_string("DERTYPE", &(Parameters.dertype),0); 
   if(errcod == IPE_KEY_NOT_FOUND) {
     Parameters.dertype = (char *) malloc(sizeof(char)*5);
     strcpy(Parameters.dertype, "NONE");
   }

   errcod = ip_string("WFN", &(Parameters.wfn),0);
   if(errcod == IPE_KEY_NOT_FOUND) {
     Parameters.wfn = (char *) malloc(sizeof(char)*5);
     strcpy(Parameters.wfn, "NONE");
   }

  if (strcmp(Parameters.dertype, "FIRST")==0 ||
      strcmp(Parameters.wfn, "DETCAS")==0) Parameters.filter_ints = 1;
  else Parameters.filter_ints = 0;

 
   /* Parameters.print_lvl is set in detci.cc */
   /* Parameters.have_special_conv is set in detci.cc */
   Parameters.ex_lvl = 2;
   Parameters.val_ex_lvl = 0;
   Parameters.maxiter = 12;
   Parameters.max_dets = 10000;
   Parameters.num_roots = 1;
   Parameters.istop = 0;
   Parameters.print_ciblks = 0;
   Parameters.S = 0;
   Parameters.opentype = PARM_OPENTYPE_UNKNOWN;
   Parameters.ref_sym = -1;
   Parameters.oei_file = PSIF_OEI;  /* always need fzc operator */
   Parameters.oei_erase = 0;
   Parameters.tei_file = PSIF_MO_TEI;
   Parameters.tei_erase = 0;
   Parameters.h0blocksize = 400;
   Parameters.h0guess_size = 100;
   Parameters.h0block_coupling = 0;
   Parameters.h0block_coupling_size = 0;
   Parameters.nprint = 20;
   Parameters.hd_ave = HD_KAVE;
   Parameters.hd_otf = TRUE;
   Parameters.nodfile = 0;
   Parameters.fzc = 1;
   Parameters.fci = 0;
   Parameters.fci_strings = 0;
   Parameters.mixed = 1;
   Parameters.mixed4 = 1;
   Parameters.r4s = 0;
   Parameters.repl_otf = 0;
   Parameters.calc_ssq = 0;
   Parameters.mpn = 0;
   Parameters.mpn_schmidt = 0;
   Parameters.wigner = 0;
   Parameters.perturbation_parameter = 1.0;
   Parameters.z_scale_H = 0;

   Parameters.ras1_lvl = -1;
   Parameters.ras1_min = -1;
   Parameters.a_ras1_lvl = -1;
   Parameters.a_ras1_min = -1;
   Parameters.b_ras1_lvl = -1;
   Parameters.b_ras1_min = -1;
   Parameters.a_ras3_max = -1;
   Parameters.b_ras3_max = -1;
   Parameters.ras3_lvl = -1;
   Parameters.ras3_max = -1;
   Parameters.ras4_lvl = -1;
   Parameters.ras4_max = -1;
   Parameters.ras34_max = -1;

   if (strcmp(Parameters.wfn, "DETCAS")==0)
     Parameters.guess_vector = PARM_GUESS_VEC_DFILE;
   else
     Parameters.guess_vector = PARM_GUESS_VEC_H0_BLOCK;

   Parameters.icore = 1;
   Parameters.diag_method = METHOD_DAVIDSON_LIU_SEM;
   Parameters.precon = PRECON_DAVIDSON;
   Parameters.update = UPDATE_DAVIDSON;
   Parameters.maxnvect = 0;
   Parameters.collapse_size = 1;
   Parameters.lse = 0;
   Parameters.lse_collapse = 3;
   Parameters.lse_tolerance = 3;
   Parameters.genci = 0;
   Parameters.neg_only = 1;
   Parameters.zero_blocks = 0;

   Parameters.nunits = 1;
   Parameters.first_tmp_unit = 50;
   Parameters.first_hd_tmp_unit = 0;
   Parameters.num_hd_tmp_units = 0;
   Parameters.first_c_tmp_unit = 0;
   Parameters.num_c_tmp_units = 0;
   Parameters.first_s_tmp_unit = 0;
   Parameters.num_s_tmp_units = 0;
   Parameters.first_d_tmp_unit = 0;
   Parameters.num_d_tmp_units = 0;
   
   Parameters.restart = 0;
   Parameters.restart_vecs = 0;
   Parameters.restart_iter = 0;
   Parameters.bendazzoli = 0;

   if (strcmp(Parameters.dertype, "FIRST")==0 ||
      strcmp(Parameters.wfn, "DETCAS")==0) {
     Parameters.convergence = 7;
     Parameters.energy_convergence = 8;
     Parameters.opdm = 1;
     Parameters.opdm_write = 1;
     Parameters.tpdm = 1;
     Parameters.tpdm_write = 1;
   }

   else {
     Parameters.convergence = 4;
     Parameters.energy_convergence = 6;
     Parameters.opdm = 0;
     Parameters.opdm_write = 0;
     Parameters.tpdm = 0;
     Parameters.tpdm_write = 0;
   }

   Parameters.opdm_file = PSIF_MO_OPDM;
   Parameters.opdm_print = 0;
   Parameters.opdm_diag = 0;
   Parameters.opdm_wrtnos = 0;
   Parameters.opdm_orbsfile = 76;
   Parameters.opdm_ave = 0;
   Parameters.opdm_orbs_root = -1;
   Parameters.tpdm_file = PSIF_MO_TPDM;
   Parameters.tpdm_print = 0;

   Parameters.nthreads = 1;
   Parameters.pthreads = 0;
   
   errcod = ip_data("EX_LVL","%d",&(Parameters.ex_lvl),0);
   errcod = ip_data("VAL_EX_LVL","%d",&(Parameters.val_ex_lvl),0);
   errcod = ip_data("MAX_DET","%d",&(Parameters.max_dets),0);
   errcod = ip_data("MAXITER","%d",&(Parameters.maxiter),0);
   errcod = ip_data("NUM_ROOTS","%d",&(Parameters.num_roots),0);
   errcod = ip_boolean("ISTOP",&(Parameters.istop),0);
   errcod = ip_data("PRINT","%d",&(Parameters.print_lvl),0);
   errcod = ip_boolean("PRINT_CIBLKS",&(Parameters.print_ciblks),0);
   errcod = ip_data("CONVERGENCE","%d",&(Parameters.convergence),0);
   errcod = ip_data("ENERGY_CONVERGENCE","%d",
               &(Parameters.energy_convergence),0);
   errcod = ip_data("S","%d",&(Parameters.S),0);

   /* this stuff was appropriate to the OPENTYPE keyword in PSI2 */
   errcod = ip_data("OPENTYPE","%s",line1,0);
   if (errcod == IPE_OK) {
      if (strcmp(line1, "NONE")==0) { 
         Parameters.opentype = PARM_OPENTYPE_NONE;
         Parameters.Ms0 = 1;
         }
      else if (strcmp(line1, "HIGHSPIN")==0) {
         Parameters.opentype = PARM_OPENTYPE_HIGHSPIN;
         Parameters.Ms0 = 0;
         }
      else if (strcmp(line1, "SINGLET")==0) {
         Parameters.opentype = PARM_OPENTYPE_SINGLET;
         Parameters.Ms0 = 1;
         }
      else Parameters.opentype = PARM_OPENTYPE_UNKNOWN;
      } 
   else { /* no opentype keyword, as appropriate for PSI3 */
     errcod = ip_data("REFERENCE","%s",line1,0);
     if (errcod == IPE_OK) {
       if (strcmp(line1, "RHF")==0) {
         Parameters.opentype = PARM_OPENTYPE_NONE;
         Parameters.Ms0 = 1;
         }
       if (strcmp(line1, "ROHF")==0) {
         if (ip_data("MULTP","%d",&tval,0) == IPE_OK) {
           if (tval == 1) {
             Parameters.opentype = PARM_OPENTYPE_SINGLET;
             Parameters.Ms0 = 1;
             }
           else {
             Parameters.opentype = PARM_OPENTYPE_HIGHSPIN;
             Parameters.Ms0 = 0;
             }
           }
         else {
           fprintf(outfile, "detci: trouble reading MULTP\n");
           exit(0);
           }
         }  /* end ROHF parsing */
       else {
         fprintf(outfile, "detci: can only handle RHF or ROHF\n");
         exit(0);
         }
       }
     else Parameters.opentype = PARM_OPENTYPE_UNKNOWN;
     } /* end PSI3 parsing */

   errcod = ip_boolean("MS0",&(Parameters.Ms0),0);
   errcod = ip_data("REF_SYM","%d",&(Parameters.ref_sym),0);
   errcod = ip_data("OEI_FILE","%d",&(Parameters.oei_file),0);
   errcod = ip_boolean("OEI_ERASE",&(Parameters.oei_erase),0);
   errcod = ip_data("TEI_FILE","%d",&(Parameters.tei_file),0);
   errcod = ip_boolean("TEI_ERASE",&(Parameters.tei_erase),0);
   errcod = ip_data("H0_BLOCKSIZE","%d",&(Parameters.h0blocksize),0);
   Parameters.h0guess_size = Parameters.h0blocksize;
   errcod = ip_data("H0_GUESS_SIZE","%d",&(Parameters.h0guess_size),0);
   if (Parameters.h0guess_size > Parameters.h0blocksize)
     Parameters.h0guess_size = Parameters.h0blocksize; 
   errcod = ip_data("H0_BLOCK_COUPLING_SIZE","%d",
                    &(Parameters.h0block_coupling_size),0);
   errcod = ip_boolean("H0_BLOCK_COUPLING",&(Parameters.h0block_coupling),0);
   errcod = ip_data("NPRINT","%d",&(Parameters.nprint),0);
   errcod = ip_boolean("FREEZE_CORE",&(Parameters.fzc),0);
   errcod = ip_boolean("FCI",&(Parameters.fci),0);
   if (Parameters.fci) Parameters.fci_strings = 1;
   errcod = ip_boolean("FCI_STRINGS",&(Parameters.fci_strings),0);
   errcod = ip_boolean("MIXED",&(Parameters.mixed),0);
   errcod = ip_boolean("MIXED4",&(Parameters.mixed4),0);
   errcod = ip_boolean("R4S",&(Parameters.r4s),0);
   errcod = ip_boolean("REPL_OTF",&(Parameters.repl_otf),0);
   errcod = ip_boolean("CALC_SSQ",&(Parameters.calc_ssq),0);
   errcod = ip_boolean("MPN",&(Parameters.mpn),0);
   if (Parameters.mpn) {
     Parameters.mpn_schmidt = FALSE;
     Parameters.wigner = TRUE;
     Parameters.guess_vector = PARM_GUESS_VEC_UNIT;
     Parameters.hd_ave = ORB_ENER;
     Parameters.update = UPDATE_DAVIDSON;
     Parameters.hd_otf = TRUE;
     Parameters.nodfile = TRUE;
     }
   errcod = ip_data("PERTURBATION_PARAMETER","%lf",
            &(Parameters.perturbation_parameter),0);
   
   if (Parameters.perturbation_parameter <= 1.0 && 
       Parameters.perturbation_parameter >= -1.0) Parameters.z_scale_H = 1;
/*
   else { fprintf(outfile, "Parameters.perturbation_parameters beyond the"
                 "bounds of -1.0 >= z <= 1.0\n");
         exit(0);
        }
*/
   errcod = ip_boolean("MPN_SCHMIDT",&(Parameters.mpn_schmidt),0);
   errcod = ip_boolean("WIGNER",&(Parameters.wigner),0);

   errcod = ip_data("A_RAS3_MAX","%d",&(Parameters.a_ras3_max),0);
   errcod = ip_data("B_RAS3_MAX","%d",&(Parameters.b_ras3_max),0);
   errcod = ip_data("RAS3_MAX","%d",&(Parameters.ras3_max),0);
   errcod = ip_data("RAS4_MAX","%d",&(Parameters.ras4_max),0);
   errcod = ip_data("RAS34_MAX","%d",&(Parameters.ras34_max),0);
   errcod = ip_data("GUESS_VECTOR", "%s", line1, 0);
   if (errcod == IPE_OK) {
      if (strcmp(line1, "UNIT")==0) 
         Parameters.guess_vector = PARM_GUESS_VEC_UNIT;
      else if (strcmp(line1, "H0_BLOCK")==0) 
         Parameters.guess_vector = PARM_GUESS_VEC_H0_BLOCK;
      else if (strcmp(line1, "DFILE")==0) 
         Parameters.guess_vector = PARM_GUESS_VEC_DFILE;
      /* else if (Parameters.mpn) Parameters.guess_vector = PARM_GUESS_VEC_UNIT; */ 
      else Parameters.guess_vector = PARM_GUESS_VEC_UNIT;
      }
   errcod = ip_data("ICORE", "%d", &(Parameters.icore),0);
   errcod = ip_boolean("GENCI",&(Parameters.genci),0) ;
   errcod = ip_data("HD_AVE","%s",line1, 0); 
   if (errcod == IPE_OK) {
     if (strcmp(line1, "HD_EXACT")==0)    Parameters.hd_ave = HD_EXACT;
     if (strcmp(line1, "HD_KAVE")==0)     Parameters.hd_ave = HD_KAVE;
     if (strcmp(line1, "ORB_ENER")==0)    Parameters.hd_ave = ORB_ENER;
     if (strcmp(line1, "EVANGELISTI")==0) Parameters.hd_ave = EVANGELISTI;
     if (strcmp(line1, "LEININGER")==0)   Parameters.hd_ave = LEININGER;
     if (strcmp(line1, "Z_KAVE")==0)      Parameters.hd_ave = Z_HD_KAVE;
     /* if (Parameters.mpn) Parameters.hd_ave = ORB_ENER; */ 
       }
   errcod = ip_boolean("HD_OTF",&(Parameters.hd_otf),0); 
   if (errcod == IPE_OK) Parameters.hd_otf = TRUE;
   errcod = ip_boolean("NODFILE",&(Parameters.nodfile),0); 
   if (Parameters.num_roots > 1) Parameters.nodfile = FALSE;

   errcod = ip_data("DIAG_METHOD","%s",line1, 0);
   if (errcod == IPE_OK) {
      if (strcmp(line1, "RSP")==0) Parameters.diag_method = METHOD_RSP;
      if (strcmp(line1, "OLSEN")==0) Parameters.diag_method = METHOD_OLSEN;
      if (strcmp(line1, "MITRUSHENKOV")==0) 
        Parameters.diag_method = METHOD_MITRUSHENKOV;
      if (strcmp(line1, "DAVIDSON")==0) 
        Parameters.diag_method = METHOD_DAVIDSON_LIU_SEM;
      if (strcmp(line1, "SEM")==0) 
        Parameters.diag_method = METHOD_DAVIDSON_LIU_SEM;
      if (strcmp(line1, "SEMTEST")==0) 
        Parameters.diag_method = METHOD_RSPTEST_OF_SEM;
      }

   errcod = ip_data("PRECONDITIONER","%s",line1, 0);
   if (errcod == IPE_OK) {
      if (strcmp(line1, "LANCZOS")==0) Parameters.precon = PRECON_LANCZOS;
      if (strcmp(line1, "DAVIDSON")==0) Parameters.precon = PRECON_DAVIDSON;
      if (strcmp(line1, "GEN_DAVIDSON")==0) 
        Parameters.precon = PRECON_GEN_DAVIDSON;
      if (strcmp(line1, "H0BLOCK")==0) Parameters.precon = PRECON_GEN_DAVIDSON;
      if (strcmp(line1, "H0BLOCK_INV")==0) 
        Parameters.precon = PRECON_H0BLOCK_INVERT;
      if (strcmp(line1, "ITER_INV")==0) 
        Parameters.precon = PRECON_H0BLOCK_ITER_INVERT;
      if (strcmp(line1, "H0BLOCK_COUPLING")==0) 
        Parameters.precon = PRECON_H0BLOCK_COUPLING;
      if (strcmp(line1, "EVANGELISTI")==0) 
        Parameters.precon = PRECON_EVANGELISTI;
      }
   errcod = ip_data("UPDATE","%s",line1, 0);
   if (errcod == IPE_OK) {
      if (strcmp(line1, "DAVIDSON")==0) Parameters.update = UPDATE_DAVIDSON;
      if (strcmp(line1, "OLSEN")==0) Parameters.update = UPDATE_OLSEN;
      }
   if (Parameters.diag_method < METHOD_DAVIDSON_LIU_SEM && 
       Parameters.update==UPDATE_DAVIDSON) {
     fprintf(outfile,"DAVIDSON update not available for OLSEN or MITRUSH"
             " iterators\n");
     Parameters.update = UPDATE_OLSEN;
     }
   if (Parameters.precon==PRECON_EVANGELISTI && (Parameters.update!=UPDATE_DAVIDSON 
        || Parameters.diag_method!=METHOD_DAVIDSON_LIU_SEM)) {
     fprintf(outfile,"EVANGELISTI preconditioner not available for OLSEN or"
                     " MITRUSH iterators or updates.\n");
     Parameters.update = UPDATE_DAVIDSON;
     }
   
   errcod = ip_boolean("ZERO_BLOCKS",&(Parameters.zero_blocks),0);
   if (Parameters.icore || !Parameters.mpn) Parameters.zero_blocks = 0;
   Parameters.num_init_vecs = Parameters.num_roots;
   errcod = ip_data("NUM_INIT_VECS","%d",&(Parameters.num_init_vecs),0);

   errcod = ip_data("COLLAPSE_SIZE", "%d", &(Parameters.collapse_size),0);
   if (Parameters.collapse_size < 1) Parameters.collapse_size = 1;

   errcod = ip_data("LSE_COLLAPSE", "%d", &(Parameters.lse_collapse),0);
   if (Parameters.lse_collapse < 1) Parameters.lse_collapse = 3;

   errcod = ip_boolean("LSE",&(Parameters.lse),0) ;

   errcod = ip_data("LSE_TOLERANCE", "%d", &(Parameters.lse_tolerance),0);

   errcod = ip_data("MAXNVECT", "%d", &(Parameters.maxnvect),0);
   if (Parameters.maxnvect == 0 &&  
       Parameters.diag_method == METHOD_DAVIDSON_LIU_SEM) {
      Parameters.maxnvect = Parameters.maxiter * Parameters.num_roots
         + Parameters.num_init_vecs;
      }
   else if (Parameters.maxnvect == 0 && 
            Parameters.diag_method == METHOD_RSPTEST_OF_SEM) {
      Parameters.maxnvect = Parameters.maxiter * Parameters.num_roots
         + Parameters.num_init_vecs;
      }
   else if (Parameters.maxnvect == 0 && 
            Parameters.diag_method == METHOD_MITRUSHENKOV) {
      Parameters.maxnvect = 2;
      }
   else if (Parameters.maxnvect == 0 && 
            Parameters.diag_method == METHOD_OLSEN) {
      Parameters.maxnvect = 1;
      }
   else { /* the user tried to specify a value for maxnvect...check it */
   /*    if (Parameters.maxnvect / (Parameters.collapse_size * 
         Parameters.num_roots) < 2) {
         fprintf(outfile, "maxnvect must be at least twice collapse_size *");
         fprintf(outfile, " num_roots.\n"); 
         exit(0);
         }
   */
      }
      
   errcod = ip_data("NUNITS", "%d", &(Parameters.nunits),0);
   errcod = ip_data("FIRST_TMP_UNIT", "%d", &(Parameters.first_tmp_unit),0);
   errcod = ip_data("FIRST_HD_TMP_UNIT","%d",&(Parameters.first_hd_tmp_unit),0);
   errcod = ip_data("FIRST_C_TMP_UNIT","%d",&(Parameters.first_c_tmp_unit),0);
   errcod = ip_data("FIRST_S_TMP_UNIT","%d",&(Parameters.first_s_tmp_unit),0);
   errcod = ip_data("FIRST_D_TMP_UNIT","%d",&(Parameters.first_d_tmp_unit),0);
   errcod = ip_data("NUM_HD_TMP_UNITS","%d",&(Parameters.num_hd_tmp_units),0);
   errcod = ip_data("NUM_C_TMP_UNITS","%d",&(Parameters.num_c_tmp_units),0);
   errcod = ip_data("NUM_S_TMP_UNITS","%d",&(Parameters.num_s_tmp_units),0);
   errcod = ip_data("NUM_D_TMP_UNITS","%d",&(Parameters.num_d_tmp_units),0);

   if (Parameters.first_hd_tmp_unit == 0) 
     Parameters.first_hd_tmp_unit = Parameters.first_tmp_unit;
/*   if ( (Parameters.num_hd_tmp_units == 0) && (!Parameters.hd_otf) ) */
   if (Parameters.num_hd_tmp_units == 0)
     Parameters.num_hd_tmp_units = 1;
   if (Parameters.first_c_tmp_unit == 0) Parameters.first_c_tmp_unit = 
      Parameters.first_hd_tmp_unit + Parameters.num_hd_tmp_units;
   if (Parameters.num_c_tmp_units == 0) Parameters.num_c_tmp_units = 
      Parameters.nunits;
   if (Parameters.first_s_tmp_unit == 0) Parameters.first_s_tmp_unit = 
      Parameters.first_c_tmp_unit + Parameters.num_c_tmp_units;
   if (Parameters.num_s_tmp_units == 0) Parameters.num_s_tmp_units = 
      Parameters.nunits;
   if (Parameters.first_d_tmp_unit == 0) Parameters.first_d_tmp_unit =
      Parameters.first_s_tmp_unit + Parameters.num_s_tmp_units;
 /*  if ( (Parameters.num_d_tmp_units == 0) && (!Parameters.nodfile) ) */
   if (Parameters.num_d_tmp_units == 0) 
     Parameters.num_d_tmp_units = 1;

   errcod = ip_boolean("RESTART",&(Parameters.restart),0);
   errcod = ip_data("RESTART_ITER","%d",&(Parameters.restart_iter),0);
   errcod = ip_data("RESTART_VECS","%d",&(Parameters.restart_vecs),0);
   if (Parameters.restart && (errcod!=IPE_OK || Parameters.restart_vecs==0)) {
      fprintf(outfile, "For RESTART must specify nonzero RESTART_VECS\n");
      exit(0);
      }
   errcod = ip_boolean("BENDAZZOLI",&(Parameters.bendazzoli),0) ;
   if (Parameters.bendazzoli & !Parameters.fci) Parameters.bendazzoli=0;

   /* Parse the OPDM stuff.  It is possible to give incompatible options,
    * but we will try to eliminate some of those.  Parameters_opdm will
    * function as the master switch for all other OPDM parameters.
    */
   errcod = ip_boolean("OPDM_PRINT",&(Parameters.opdm_print),0);
   errcod = ip_data("OPDM_FILE","%d",&(Parameters.opdm_file),0);
   errcod = ip_boolean("WRTNOS",&(Parameters.opdm_wrtnos),0);
   errcod = ip_boolean("OPDM_DIAG",&(Parameters.opdm_diag),0);
   errcod = ip_boolean("OPDM_AVE",&(Parameters.opdm_ave),0);
   errcod = ip_data("ORBSFILE","%d",&(Parameters.opdm_orbsfile),0);
   errcod = ip_data("ORBS_ROOT","%d",&(Parameters.opdm_orbs_root),0);
   
   if (Parameters.opdm_orbs_root != -1) Parameters.opdm_orbs_root -= 1;
   if (Parameters.opdm_orbs_root < 0) Parameters.opdm_orbs_root = 0;
   if (Parameters.opdm_wrtnos) Parameters.opdm_diag = 1;
   if (Parameters.opdm_print || Parameters.opdm_diag || Parameters.opdm_wrtnos 
       || Parameters.opdm_ave) Parameters.opdm = 1;
   errcod = ip_boolean("OPDM",&(Parameters.opdm),0);
   if (Parameters.opdm) Parameters.opdm_write = 1;
   errcod = ip_boolean("OPDM_WRITE",&(Parameters.opdm_write),0);
   errcod = ip_boolean("OPDM_PRINT",&(Parameters.opdm_print),0);
   errcod = ip_data("OPDM_FILE","%d",&(Parameters.opdm_file),0);
   errcod = ip_data("OPDM_DIAG","%d",&(Parameters.opdm_diag),0);
   errcod = ip_data("WRTNOS","%d",&(Parameters.opdm_wrtnos),0);
   errcod = ip_data("OPDM_AVE", "%d", &(Parameters.opdm_ave),0);
   errcod = ip_data("ORBSFILE", "%d", &(Parameters.opdm_orbsfile),0);
   errcod = ip_data("ORBS_ROOT", "%d", &(Parameters.opdm_orbs_root),0);
   if (Parameters.opdm_orbs_root == -1) 
     Parameters.opdm_orbs_root = 0;
    /* Parameters.opdm_orbs_root = Parameters.num_roots-1; */
   else Parameters.opdm_orbs_root -= 1;

   errcod = ip_boolean("TPDM",&(Parameters.tpdm),0);
   if (Parameters.tpdm) Parameters.tpdm_write = 1;
   errcod = ip_boolean("TPDM_WRITE",&(Parameters.tpdm_write),0);
   errcod = ip_boolean("TPDM_PRINT",&(Parameters.tpdm_print),0);
   errcod = ip_data("TPDM_FILE","%d",&(Parameters.tpdm_file),0);

   if (Parameters.guess_vector == PARM_GUESS_VEC_DFILE &&
       strcmp(Parameters.wfn, "DETCAS")!=0) {
      file30_init();
      i = file30_rd_phase_check();
      file30_close();
      if (!i) {
         fprintf(outfile, "Can't use d file guess: SCF phase not checked\n");
         if (Parameters.h0guess_size) {
            Parameters.guess_vector = PARM_GUESS_VEC_H0_BLOCK;
            if (Parameters.precon == PRECON_GEN_DAVIDSON)
              Parameters.precon = PRECON_H0BLOCK_ITER_INVERT;
         }
         else Parameters.guess_vector = PARM_GUESS_VEC_UNIT;
      }
     }
   if (Parameters.num_init_vecs < Parameters.num_roots)
     Parameters.num_init_vecs = Parameters.num_roots;
   if (Parameters.guess_vector == PARM_GUESS_VEC_UNIT &&
       Parameters.num_init_vecs > 1) {
     Parameters.guess_vector = PARM_GUESS_VEC_H0_BLOCK;
     fprintf(outfile,"Warning: Unit vec option not available for more than"
             " one root\n");
     }
   if (Parameters.guess_vector == PARM_GUESS_VEC_UNIT)
     Parameters.h0blocksize = Parameters.h0guess_size = 1;

   errcod = ip_data("NTHREADS", "%d", &(Parameters.nthreads),0);
   if (Parameters.nthreads < 1) Parameters.nthreads = 1;
   errcod = ip_boolean("PTHREADS",&(Parameters.pthreads),0);
   if (!Parameters.pthreads) Parameters.nthreads = 1;
}


/*
** print_parameters(): Function prints the program's running parameters
**   found in the Parameters structure.
*/
void print_parameters(void)
{
   fprintf(outfile, "\n");
   fprintf(outfile, "PARAMETERS: \n");
   fprintf(outfile, "   EX LVL        =   %6d      H0 BLOCKSIZE =   %6d\n", 
      Parameters.ex_lvl, Parameters.h0blocksize);
   fprintf(outfile, "   VAL EX LVL    =   %6d      H0 GUESS SIZE=   %6d\n", 
      Parameters.val_ex_lvl, Parameters.h0guess_size);
   fprintf(outfile, "   H0COUPLINGSIZE=   %6d      H0 COUPLING  =   %6s\n", 
      Parameters.h0block_coupling_size, Parameters.h0block_coupling ? "yes" : "no");
   fprintf(outfile, "   NPRINT        =   %6d      MAX DET      =   %6d\n", 
      Parameters.nprint, Parameters.max_dets);
   fprintf(outfile, "   MAXITER       =   %6d      FREEZE CORE  =   %6s\n", 
      Parameters.maxiter, Parameters.fzc ? "yes" : "no");
   fprintf(outfile, "   NUM ROOTS     =   %6d      ICORE        =   %6d\n", 
      Parameters.num_roots, Parameters.icore);
   fprintf(outfile, "   PRINT         =   %6d      FCI          =   %6s\n", 
      Parameters.print_lvl, Parameters.fci ? "yes" : "no");
   if (Parameters.have_special_conv) 
      fprintf(outfile, 
         "   CONV          =   %8.2g    MIXED        =   %6s\n", 
         Parameters.special_conv, Parameters.mixed ? "yes" : "no");
   else
      fprintf(outfile, "   CONV          =   %6d      MIXED        =   %6s\n", 
         Parameters.convergence, Parameters.mixed ? "yes" : "no");

   fprintf(outfile, "   E CONV        =   %6d      MIXED4       =   %6s\n", 
      Parameters.energy_convergence, Parameters.mixed4 ? "yes" : "no");
   fprintf(outfile, "   OEI FILE      =   %6d      R4S          =   %6s\n", 
      Parameters.oei_file, Parameters.r4s ? "yes" : "no");
   fprintf(outfile, "   OEI ERASE     =   %6s      REPL OTF     =   %6s\n",  
      Parameters.oei_erase ? "yes" : "no", Parameters.repl_otf ? "yes" : "no");
   fprintf(outfile, "   TEI FILE      =   %6d      DIAG METHOD  =   ", 
      Parameters.tei_file);

   switch (Parameters.diag_method) {
      case 0:
         fprintf(outfile, "%6s\n", "RSP");
         break;
      case 1:
         fprintf(outfile, "%6s\n", "OLSEN");
         break;
      case 2:
         fprintf(outfile, "%6s\n", "MITRUS");
         break;
      case 3:
         fprintf(outfile, "%6s\n", "SEM");
         break;
      case 4:
         fprintf(outfile, "%6s\n", "SEMTEST");
         break;
      default:
         fprintf(outfile, "%6s\n", "???");
         break;
      } 

   fprintf(outfile, "   PRECONDITIONER= ");
   switch (Parameters.precon) {
      case PRECON_LANCZOS:
         fprintf(outfile, "%6s", " LANCZOS    ");
         break;
      case PRECON_DAVIDSON:
         fprintf(outfile, "%6s", "DAVIDSON    ");
         break;
      case PRECON_GEN_DAVIDSON:
         fprintf(outfile, "%6s", "GEN_DAVIDSON");
         break;
      case PRECON_H0BLOCK_INVERT:
         fprintf(outfile, "%6s", "H0BLOCK_INV ");
         break;
      case PRECON_H0BLOCK_ITER_INVERT:
         fprintf(outfile, "%6s", "ITER_INV    ");
         break;
      case PRECON_H0BLOCK_COUPLING:
         fprintf(outfile, "%6s", "H0_COUPLING ");
         break;
      case PRECON_EVANGELISTI:
         fprintf(outfile, "%6s", "EVANGELISTI ");
         break;
      default:
         fprintf(outfile, "%6s", "???         ");
         break;
      } 

   fprintf(outfile, "  UPDATE       =   ");
   switch (Parameters.update) {
     case 1:
       fprintf(outfile, "%6s\n", "DAVIDSON");
       break;
     case 2:
       fprintf(outfile, "%6s\n", "OLSEN");  
       break;
     default:
       fprintf(outfile, "%6s\n", "???");
       break;
      }

   fprintf(outfile, "   S             =   %6d      Ms0          =   %6s\n",
      Parameters.S, Parameters.Ms0 ? "yes" : "no");           
   fprintf(outfile, "   TEI ERASE     =   %6s      MAXNVECT     =   %6d\n", 
      Parameters.tei_erase ? "yes" : "no", Parameters.maxnvect);
   fprintf(outfile, "   RESTART       =   %6s      RESTART VECS =   %6d\n",
      Parameters.restart ? "yes" : "no", Parameters.restart_vecs);
   fprintf(outfile, "   GUESS VECTOR  =  ");
   switch (Parameters.guess_vector) {
      case PARM_GUESS_VEC_UNIT:
         fprintf(outfile, "%7s", "UNIT");
         break; 
      case PARM_GUESS_VEC_H0_BLOCK:
         fprintf(outfile, "%7s", "H0BLOCK");
         break;
      case PARM_GUESS_VEC_DFILE:
         fprintf(outfile, "%7s", "D FILE");
         break;
      default:
         fprintf(outfile, "%7s", "???");
         break;
      }
   fprintf(outfile, "      OPENTYPE     = ");
   switch (Parameters.opentype) {
      case PARM_OPENTYPE_NONE:
         fprintf(outfile, "%8s\n", "NONE");
         break;
      case PARM_OPENTYPE_HIGHSPIN:
         fprintf(outfile, "%8s\n", "HIGHSPIN");
         break;
      case PARM_OPENTYPE_SINGLET:
         fprintf(outfile, "%8s\n", "SINGLET");
         break;
      default:
         fprintf(outfile, "%8s\n", "???"); 
         break;
      }
   fprintf(outfile, "   GENCI ALG     =   %6s",
      Parameters.genci ? "yes" : "no");
   if (Parameters.ref_sym == -1)
      fprintf(outfile, "      REF SYM      =   %6s\n", "auto");
   else
      fprintf(outfile, "      REF SYM      =   %6d\n", Parameters.ref_sym);

   fprintf(outfile, "   COLLAPSE SIZE =   %6d", Parameters.collapse_size);
   fprintf(outfile, "      HD AVE       =");
   switch (Parameters.hd_ave) {
     case HD_EXACT:
       fprintf(outfile," %11s\n", "HD_EXACT");
       break;
     case HD_KAVE:
       fprintf(outfile," %11s\n", "HD_KAVE");
       break;
     case ORB_ENER:
       fprintf(outfile," %11s\n", "ORB_ENER");
       break;
     case EVANGELISTI:
       fprintf(outfile," %11s\n", "EVANGELISTI");
       break;
     case LEININGER:
       fprintf(outfile," %11s\n", "LEININGER");
       break;
     default:
       fprintf(outfile," %11s\n", "???");       
       break;
     }

   fprintf(outfile, "   LSE           =   %6s      LSE ITER     =   %6d\n", 
           Parameters.lse ? "yes" : "no", Parameters.lse_iter);
   fprintf(outfile, "   HD OTF        =   %6s      NO DFILE     =   %6s\n", 
           Parameters.hd_otf ? "yes" : "no", Parameters.nodfile ? "yes":"no");
   fprintf(outfile, "   MPN           =   %6s      MPN SCHMIDT  =   %6s\n",
           Parameters.mpn ? "yes":"no", Parameters.mpn_schmidt ? "yes":"no");
   fprintf(outfile, "   WIGNER        =   %6s      ZERO BLOCKS  =   %6s\n", 
           Parameters.wigner ? "yes":"no", Parameters.zero_blocks ? "yes":"no");
   fprintf(outfile, "   PERT Z        =   %1.4f      NTHREADS     =        %d\n",
           Parameters.perturbation_parameter, Parameters.nthreads);
   fprintf(outfile, "   PTHREADS      =   %6s\n",
           Parameters.pthreads ? "yes":"no");
   fprintf(outfile, "\n   FILES         =     %3d %3d %3d %3d\n",
      Parameters.first_hd_tmp_unit, Parameters.first_c_tmp_unit,
      Parameters.first_s_tmp_unit, Parameters.first_d_tmp_unit);

   fprintf(outfile, "\n") ;
   fflush(outfile) ;
}


/*
** set_ras_parms(): Set the RAS parameters or their conventional equivalents
**   (i.e. fermi level, etc).
**
*/
void set_ras_parms(void)
{
   int i,j,cnt;
   int errcod;
   int tot_expl_el,nras2alp,nras2bet,betsocc;
   int *ras1, *ras2, *ras3;
   int *orbsym;

   for (i=0,j=0; i<CalcInfo.nirreps; i++) j += CalcInfo.ras_opi[0][i];
   Parameters.a_ras1_lvl = Parameters.b_ras1_lvl = Parameters.ras1_lvl = j-1;

   /* figure out how many electrons are in RAS II */
   /* alpha electrons */
   for (i=0,nras2alp=0,betsocc=0; i<CalcInfo.nirreps; i++) {
      j = CalcInfo.docc[i] - CalcInfo.frozen_docc[i] - CalcInfo.ras_opi[0][i];
      if (Parameters.opentype == PARM_OPENTYPE_HIGHSPIN) {
         j += CalcInfo.socc[i];
         }
      else if (Parameters.opentype == PARM_OPENTYPE_SINGLET) {
         if (betsocc + CalcInfo.socc[i] <= CalcInfo.spab)
            betsocc += CalcInfo.socc[i];
         else {
            j += CalcInfo.socc[i] - (CalcInfo.spab - betsocc);
            betsocc = CalcInfo.spab;
            }
         }
      if (j > 0) nras2alp += j;
      if (j > CalcInfo.ras_opi[1][i]) {
         fprintf(outfile, "(set_ras_parms): detecting %d electrons ",
            j - CalcInfo.ras_opi[1][i]);
         fprintf(outfile, "in RAS III for irrep %d.\n", i);
         fprintf(outfile, "Some parts of DETCI assume all elec in I and II\n");
         }
      }
   /* beta electrons */
   for (i=0,nras2bet=0,betsocc=0; i<CalcInfo.nirreps; i++) {
      j = CalcInfo.docc[i] - CalcInfo.frozen_docc[i] - CalcInfo.ras_opi[0][i];
      if (Parameters.opentype == PARM_OPENTYPE_SINGLET && CalcInfo.socc[i]) {
         if (betsocc + CalcInfo.socc[i] <= CalcInfo.spab)
            j += CalcInfo.socc[i];
         else {
            j += CalcInfo.spab - betsocc;
            betsocc = CalcInfo.spab;
            }
         }
      if (j > 0) nras2bet += j;
      if (j > CalcInfo.ras_opi[1][i]) {
         fprintf(outfile, "(set_ras_parms): detecting %d electrons ",
            j - CalcInfo.ras_opi[1][i]);
         fprintf(outfile, "in RAS III for irrep %d.\n", i);
         fprintf(outfile, "Some parts of DETCI assume all elec in I and II\n");
         }
      }

   Parameters.a_ras1_max = (CalcInfo.num_alp_expl >
         Parameters.a_ras1_lvl + 1) ? Parameters.a_ras1_lvl + 1 :
         (CalcInfo.num_alp_expl) ;
   if (Parameters.fzc) Parameters.a_ras1_max += CalcInfo.num_fzc_orbs;

   Parameters.b_ras1_max = (CalcInfo.num_bet_expl >
         Parameters.b_ras1_lvl + 1) ? Parameters.b_ras1_lvl + 1:
         (CalcInfo.num_bet_expl) ;
   if (Parameters.fzc) Parameters.b_ras1_max += CalcInfo.num_fzc_orbs;

   for (i=0,j=0; i<CalcInfo.nirreps; i++) j += CalcInfo.ras_opi[1][i];
   Parameters.ras3_lvl = Parameters.ras1_lvl + j + 1;

   for (i=0,j=0; i<CalcInfo.nirreps; i++) j += CalcInfo.ras_opi[2][i];
   Parameters.ras4_lvl = Parameters.ras3_lvl + j;


   /* check Parameters to make sure everything consistent */

   /* deduce Parameters.a_ras3_max and Parameters.b_ras3_max if needed */
   if (Parameters.a_ras3_max == -1 || Parameters.b_ras3_max == -1) {
      if (Parameters.ras3_max != -1) { /* have parsed ras3_max */
         Parameters.a_ras3_max = (Parameters.ras3_max <= CalcInfo.num_alp_expl)
            ? Parameters.ras3_max : CalcInfo.num_alp_expl;
         Parameters.b_ras3_max = (Parameters.ras3_max <= CalcInfo.num_bet_expl)
            ? Parameters.ras3_max : CalcInfo.num_bet_expl;
         }
      else {
         Parameters.a_ras3_max = (Parameters.ex_lvl <= CalcInfo.num_alp_expl) 
            ? Parameters.ex_lvl : CalcInfo.num_alp_expl;
         Parameters.b_ras3_max = (Parameters.ex_lvl <= CalcInfo.num_bet_expl) 
            ? Parameters.ex_lvl : CalcInfo.num_bet_expl; 
         }
      }

   if (Parameters.ras4_max != -1) { /* have parsed */
      Parameters.a_ras4_max = (Parameters.ras4_max <= CalcInfo.num_alp_expl)
        ? Parameters.ras4_max : CalcInfo.num_alp_expl;
      Parameters.b_ras4_max = (Parameters.ras4_max <= CalcInfo.num_bet_expl)
        ? Parameters.ras4_max : CalcInfo.num_bet_expl;
      }
   else {
      Parameters.a_ras4_max = Parameters.a_ras3_max;
      Parameters.b_ras4_max = Parameters.b_ras3_max;
   }

   if (Parameters.ras34_max != -1) { /* have parsed */
      Parameters.a_ras34_max = Parameters.ras34_max;
      Parameters.b_ras34_max = Parameters.ras34_max;
      }
   else {
      Parameters.a_ras34_max = Parameters.a_ras3_max;
      Parameters.b_ras34_max = Parameters.b_ras3_max;
      }

   i = Parameters.ras4_lvl - Parameters.ras3_lvl;
   if (Parameters.a_ras3_max > i) Parameters.a_ras3_max = i;
   if (Parameters.b_ras3_max > i) Parameters.b_ras3_max = i;

   i = CalcInfo.num_ci_orbs - Parameters.ras4_lvl;
   if (Parameters.a_ras4_max > i) Parameters.a_ras4_max = i;
   if (Parameters.b_ras4_max > i) Parameters.b_ras4_max = i;

   i = CalcInfo.num_ci_orbs - Parameters.ras3_lvl;
   if (Parameters.a_ras34_max > i) Parameters.a_ras34_max = i;
   if (Parameters.b_ras34_max > i) Parameters.b_ras34_max = i;

   i = (CalcInfo.num_alp_expl <= Parameters.a_ras1_lvl + 1) ? 
      CalcInfo.num_alp_expl : Parameters.a_ras1_lvl + 1;
   Parameters.a_ras1_min = i - Parameters.ex_lvl -
      Parameters.val_ex_lvl;
   if (Parameters.a_ras1_min < 0) Parameters.a_ras1_min = 0;
   Parameters.a_ras1_min += CalcInfo.num_fzc_orbs;
   Parameters.a_ras1_min += CalcInfo.num_cor_orbs;

   i = (CalcInfo.num_bet_expl <= Parameters.b_ras1_lvl + 1) ? 
      CalcInfo.num_bet_expl : Parameters.b_ras1_lvl + 1;
   Parameters.b_ras1_min = i - Parameters.ex_lvl -
      Parameters.val_ex_lvl;
   if (Parameters.b_ras1_min < 0) Parameters.b_ras1_min = 0;
   Parameters.b_ras1_min += CalcInfo.num_fzc_orbs;
   Parameters.b_ras1_min += CalcInfo.num_cor_orbs;

   tot_expl_el = CalcInfo.num_alp_expl + CalcInfo.num_bet_expl;
   if (Parameters.ras3_max == -1) {
      Parameters.ras3_max = (Parameters.ex_lvl <= tot_expl_el) ?
         Parameters.ex_lvl : tot_expl_el ;
      }
   else {
      if (Parameters.ras3_max > tot_expl_el) 
         Parameters.ras3_max = tot_expl_el;
      }

   if (Parameters.ras34_max == -1) Parameters.ras34_max = Parameters.ras3_max;

   i = 2 * (Parameters.ras4_lvl - Parameters.ras3_lvl);
   if (i < Parameters.ras3_max) Parameters.ras3_max = i;

   i = 2 * (CalcInfo.num_ci_orbs - Parameters.ras3_lvl);
   if (i < Parameters.ras34_max) Parameters.ras34_max = i;

   i = (tot_expl_el < 2*(Parameters.ras1_lvl + 1)) ? tot_expl_el :
      2*(Parameters.ras1_lvl + 1) ;
   
   Parameters.ras1_min = i - Parameters.ex_lvl - 
      Parameters.val_ex_lvl + 2 * CalcInfo.num_fzc_orbs;
   
   if (Parameters.a_ras1_min + Parameters.b_ras1_min > Parameters.ras1_min)
      Parameters.ras1_min = Parameters.a_ras1_min + Parameters.b_ras1_min;

   if (Parameters.ras4_max == -1) {
      Parameters.ras4_max = (Parameters.ex_lvl <= tot_expl_el) ?
         Parameters.ex_lvl : tot_expl_el;
      }

   i = 2 * (CalcInfo.num_ci_orbs - Parameters.ras4_lvl);
   if (i < Parameters.ras4_max) Parameters.ras4_max = i;

   if (Parameters.print_lvl) {
      fprintf(outfile, "   RAS1 LVL     =   %6d      A RAS3 MAX   =   %6d\n",
         Parameters.ras1_lvl, Parameters.a_ras3_max);
      fprintf(outfile, "   RAS1 MIN     =   %6d      B RAS3 MAX   =   %6d\n",
         Parameters.ras1_min, Parameters.b_ras3_max);
      fprintf(outfile, "   A RAS1 LVL   =   %6d      RAS4 LVL     =   %6d\n", 
         Parameters.a_ras1_lvl, Parameters.ras4_lvl);
      fprintf(outfile, "   A RAS1 MIN   =   %6d      A RAS4 MAX   =   %6d\n", 
         Parameters.a_ras1_min, Parameters.a_ras4_max);
      fprintf(outfile, "   A RAS1 MAX   =   %6d      B RAS4 MAX   =   %6d\n", 
         Parameters.a_ras1_max, Parameters.b_ras4_max);
      fprintf(outfile, "   B RAS1 LVL   =   %6d      RAS4 MAX     =   %6d\n", 
         Parameters.b_ras1_lvl, Parameters.ras4_max);
      fprintf(outfile, "   B RAS1 MIN   =   %6d      A RAS34 MAX  =   %6d\n", 
         Parameters.b_ras1_min, Parameters.a_ras34_max);
      fprintf(outfile, "   B RAS1 MAX   =   %6d      B RAS34 MAX  =   %6d\n", 
         Parameters.b_ras1_max, Parameters.b_ras34_max);
      fprintf(outfile, "   RAS3 LVL     =   %6d      RAS34 MAX    =   %6d\n", 
         Parameters.ras3_lvl, Parameters.ras34_max);
      fprintf(outfile, "   RAS3 MAX     =   %6d\n", Parameters.ras3_max);

   
      fprintf(outfile, "\n");
      fprintf(outfile, "   DOCC         = ") ;
      for (i=0; i<CalcInfo.nirreps; i++) {
         fprintf(outfile, "%2d ", CalcInfo.docc[i]) ;
         }
      fprintf(outfile, "\n   SOCC         = ") ;
      for (i=0; i<CalcInfo.nirreps; i++) {
         fprintf(outfile, "%2d ", CalcInfo.socc[i]) ;
         }
      fprintf(outfile, "\n   FROZEN DOCC  = ") ;
      for (i=0; i<CalcInfo.nirreps; i++) {
         fprintf(outfile, "%2d ", CalcInfo.frozen_docc[i]) ;
         }
      fprintf(outfile, "\n   FROZEN UOCC  = ") ;
      for (i=0; i<CalcInfo.nirreps; i++) {
         fprintf(outfile, "%2d ", CalcInfo.frozen_uocc[i]) ;
         }
      fprintf(outfile, "\n");
      for (i=0; i<4; i++) {
         fprintf(outfile, "   RAS %d        = ",i+1);
         for (j=0; j<CalcInfo.nirreps; j++) {
            fprintf(outfile,"%2d ",CalcInfo.ras_opi[i][j]);
            }
         fprintf(outfile, "\n");
         }

      fprintf(outfile,
         "*******************************************************\n\n");
      }
}



