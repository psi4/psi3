#include <stdio.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <file30.h>
#include <qt.h>
#include "structs.h"
#define EXTERN
#include "globals.h"

/*
** GET_MO_INFO
** 
** Reads file30 & input.dat and gets all sorts of useful information about 
** the molecular orbitals (such as their reordering array, the docc
** array, frozen orbitals, etc.)
**
** Created by C. David Sherrill on 17 November 1994
**
** Updated
** CDS  1/18/95 to read SCF eigenvalues also (for MP2 guess vector)
** CDS  1/ 5/97 to use nifty new ras_set() function (which transqt has been
**              using for some time).
**
*/
void get_mo_info(void)
{
   int i, j, k, tmp, cnt, irrep, errcod, errbad;
   int size;
   double *eig_unsrt;
   int parsed_ras1=0, parsed_ras2=0, do_ras4;

   CalcInfo.maxKlist = 0.0; 
   
   file30_init();
   CalcInfo.nirreps = file30_rd_nirreps();
   CalcInfo.nso = file30_rd_nmo();
   CalcInfo.nmo = file30_rd_nmo();
   CalcInfo.iopen = file30_rd_iopen();
   CalcInfo.labels = file30_rd_irr_labs();
   CalcInfo.orbs_per_irr = file30_rd_orbspi();
   CalcInfo.so_per_irr = file30_rd_sopi();
   CalcInfo.closed_per_irr = file30_rd_clsdpi();
   CalcInfo.open_per_irr = file30_rd_openpi();
   CalcInfo.enuc = file30_rd_enuc();
   CalcInfo.escf = file30_rd_escf();
   CalcInfo.efzc = file30_rd_efzc();
   eig_unsrt = file30_rd_evals();
   file30_close();
 
   if (CalcInfo.iopen && Parameters.opentype == PARM_OPENTYPE_NONE) {
      fprintf(outfile, "Warning: iopen=1,opentype=none. Making iopen=0\n");
      CalcInfo.iopen = 0;
      }
   else if (!CalcInfo.iopen && (Parameters.opentype == PARM_OPENTYPE_HIGHSPIN
      || Parameters.opentype == PARM_OPENTYPE_SINGLET)) {
      fprintf(outfile,"Warning: iopen=0,opentype!=closed. Making iopen=1\n");
      CalcInfo.iopen = 1;
      }
   if (Parameters.ref_sym >= CalcInfo.nirreps) {
      fprintf(outfile,"Warning: ref_sym >= nirreps.  Setting ref_sym=0\n");
      Parameters.ref_sym = 0;
      }

   CalcInfo.docc = init_int_array(CalcInfo.nirreps);
   CalcInfo.socc = init_int_array(CalcInfo.nirreps);
   CalcInfo.frozen_docc = init_int_array(CalcInfo.nirreps);
   CalcInfo.frozen_uocc = init_int_array(CalcInfo.nirreps);
   CalcInfo.reorder = init_int_array(CalcInfo.nmo);
   CalcInfo.ras_opi = init_int_matrix(4,CalcInfo.nirreps);
      
   if (!ras_set(CalcInfo.nirreps, CalcInfo.nmo, Parameters.fzc, 
                CalcInfo.orbs_per_irr, CalcInfo.docc, CalcInfo.socc, 
                CalcInfo.frozen_docc, CalcInfo.frozen_uocc, 
                CalcInfo.ras_opi, CalcInfo.reorder, 1)) 
   { 
     fprintf(outfile, "Error in ras_set().  Aborting.\n");
     exit(1);
   }
   
   /* calculate number of orbitals active in CI */
   CalcInfo.num_ci_orbs = CalcInfo.nmo ;
   for (i=0; i<CalcInfo.nirreps; i++) {
      CalcInfo.num_ci_orbs -= CalcInfo.frozen_uocc[i] ;
      }

   if ((CalcInfo.num_ci_orbs * (CalcInfo.num_ci_orbs + 1)) / 2 > IOFF_MAX) {
      fprintf(outfile, "Error: IOFF_MAX not large enough!\n");
      exit(1);
   }

   /* Compute maximum number of orbitals per irrep including
   ** and not including fzv
   */
  CalcInfo.max_orbs_per_irrep = 0;
  CalcInfo.max_pop_per_irrep = 0;
  for (i=0; i<CalcInfo.nirreps; i++) {
     if (CalcInfo.max_orbs_per_irrep < CalcInfo.orbs_per_irr[i])
       CalcInfo.max_orbs_per_irrep = CalcInfo.orbs_per_irr[i];
     if (CalcInfo.max_pop_per_irrep < (CalcInfo.orbs_per_irr[i] - 
                                   CalcInfo.frozen_uocc[i]))
       CalcInfo.max_pop_per_irrep = CalcInfo.orbs_per_irr[i] -
                                    CalcInfo.frozen_uocc[i];      
     }


   /* construct the "ordering" array, which maps the other direction */
   /* i.e. from a CI orbital to a Pitzer orbital                     */
   CalcInfo.order = init_int_array(CalcInfo.nmo);
   for (i=0; i<CalcInfo.nmo; i++) {
      j = CalcInfo.reorder[i];
      CalcInfo.order[j] = i;
      }


   if (Parameters.print_lvl > 4) {
      fprintf(outfile, "\nReordering array = \n");
      for (i=0; i<CalcInfo.nmo; i++) {
         fprintf(outfile, "%3d ", CalcInfo.reorder[i]);
         }
      fprintf(outfile, "\n");
      }

   CalcInfo.nmotri = (CalcInfo.nmo * (CalcInfo.nmo + 1)) / 2 ;

   /* transform orbsym vector to new MO order */
   CalcInfo.orbsym = init_int_array(CalcInfo.nmo);
   CalcInfo.scfeigval = init_array(CalcInfo.nmo);

   for (i=0,cnt=0; i<CalcInfo.nirreps; i++) {
      for (j=0; j<CalcInfo.orbs_per_irr[i]; j++,cnt++) {
         k = CalcInfo.reorder[cnt];
         CalcInfo.orbsym[k] = i;
         }
      }

   for (i=0; i<CalcInfo.nmo; i++) {
      j = CalcInfo.reorder[i];
      CalcInfo.scfeigval[j] = eig_unsrt[i];
      }
   free(eig_unsrt);

   /* calculate number of electrons */
   CalcInfo.num_alp = CalcInfo.num_bet = CalcInfo.spab = 0;
   if (Parameters.opentype == PARM_OPENTYPE_NONE ||
       Parameters.opentype == PARM_OPENTYPE_HIGHSPIN) {
      for (i=0; i<CalcInfo.nirreps; i++) {
         CalcInfo.num_alp += CalcInfo.docc[i] + CalcInfo.socc[i];
         CalcInfo.num_bet += CalcInfo.docc[i];
         }
      }
   else if (Parameters.opentype == PARM_OPENTYPE_SINGLET) {
      for (i=0; i<CalcInfo.nirreps; i++) { /* closed-shell part */
         CalcInfo.spab += CalcInfo.socc[i];      
         CalcInfo.num_alp += CalcInfo.docc[i];
         CalcInfo.num_bet += CalcInfo.docc[i];
         }
      if (CalcInfo.spab % 2) { 
         fprintf(outfile,"For opentype=singlet must have even number ");
         fprintf(outfile,"of socc electrons!\n");
         exit(1);
         }
      CalcInfo.spab /= 2;
      tmp = 0;
      for (i=0; i<CalcInfo.nirreps; i++) {
         j = CalcInfo.socc[i];
         k = 0;
         while (k < j) {
            if (tmp < CalcInfo.spab) {
               CalcInfo.num_alp++;
               tmp++; 
               k++;
               }
            else { 
               CalcInfo.num_bet++;
               tmp++;
               k++;
               } 
            }
         }
      }
   else {
      fprintf(outfile, "(get_mo_info): Can't handle opentype = %d\n",
         Parameters.opentype);
      exit(1);
      }

   CalcInfo.num_fzv_orbs = 0;
   for (i=0; i<CalcInfo.nirreps; i++) 
      CalcInfo.num_fzv_orbs += CalcInfo.frozen_uocc[i];  

   CalcInfo.num_fzc_orbs = 0;
   CalcInfo.num_cor_orbs = 0;
   if (Parameters.fzc) {
      for (i=0; i<CalcInfo.nirreps; i++) {
         j = CalcInfo.frozen_docc[i];
         CalcInfo.num_ci_orbs -= j;
         CalcInfo.num_fzc_orbs += j;
         }
      }
   else {
      for (i=0; i<CalcInfo.nirreps; i++) {
         CalcInfo.num_cor_orbs += CalcInfo.frozen_docc[i];
         } 
      }

   CalcInfo.num_alp_expl = CalcInfo.num_alp - CalcInfo.num_fzc_orbs;
   CalcInfo.num_bet_expl = CalcInfo.num_bet - CalcInfo.num_fzc_orbs;


   /* construct the CalcInfo.ras_orbs array (may not be of any use now) */
   cnt = 0;
   for (i=0; i<4; i++) {
     CalcInfo.ras_orbs[i] = init_int_matrix(CalcInfo.nirreps,
       CalcInfo.num_ci_orbs);
     for (irrep=0; irrep<CalcInfo.nirreps; irrep++) {
       for (j=0; j<CalcInfo.ras_opi[i][irrep]; j++) {
         CalcInfo.ras_orbs[i][irrep][j] = cnt++;
       }
     }
   }

   if (Parameters.print_lvl > 0) {
      fprintf(outfile, "ORBITALS:\n") ;
      /*
      fprintf(outfile, "   DOCC         = ") ;
      for (i=0; i<CalcInfo.nirreps; i++) {
         fprintf(outfile, "%d ", CalcInfo.docc[i]) ;
         }
      fprintf(outfile, "\n   SOCC         = ") ;
      for (i=0; i<CalcInfo.nirreps; i++) {
         fprintf(outfile, "%d ", CalcInfo.socc[i]) ;
         }
      fprintf(outfile, "\n   FROZEN_DOCC  = ") ;
      for (i=0; i<CalcInfo.nirreps; i++) {
         fprintf(outfile, "%d ", CalcInfo.frozen_docc[i]) ;
         }
      fprintf(outfile, "\n   FROZEN_UOCC  = ") ;
      for (i=0; i<CalcInfo.nirreps; i++) {
         fprintf(outfile, "%d ", CalcInfo.frozen_uocc[i]) ;
         }
      fprintf(outfile, "\n");
      */
      fprintf(outfile, "   NMO          =   %6d      NUM ALP      =   %6d\n",
         CalcInfo.nmo, CalcInfo.num_alp);
      fprintf(outfile, "   ORBS IN CI   =   %6d      NUM ALP EXPL =   %6d\n",
         CalcInfo.num_ci_orbs, CalcInfo.num_alp_expl);
      fprintf(outfile, "   FROZEN CORE  =   %6d      NUM BET      =   %6d\n",
         CalcInfo.num_fzc_orbs, CalcInfo.num_bet);
      fprintf(outfile, "   RESTR CORE   =   %6d      NUM BET EXPL =   %6d\n",
         CalcInfo.num_cor_orbs, CalcInfo.num_bet_expl);
      fprintf(outfile, "   IOPEN        =   %6s\n", CalcInfo.iopen ? "yes" :
         "no");
      }
}

