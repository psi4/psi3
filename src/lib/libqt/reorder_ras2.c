#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>

/*
** reorder_ras2()
**
** This function constructs a reordering array appropriate for RAS 
** wavefunctions.  The reordering array takes a basis function in 
** Pitzer ordering (orbitals grouped according to irrep) and gives the
** corresponding index in the RAS numbering scheme.
**
** C. David Sherrill
** Center for Computational Quantum Chemistry
** University of Georgia, 1995
**
** Updated 8/16/95 by CDS
**    If there is a RAS3 vector, then assume there is a RAS IV space (should
**    be useful for limiting CISD[TQ]'s, etc).
** Updated 4/24/96 by CDS
**    Use a new (simpler) ordering scheme to facilitate on-the-fly addressing
**    of strings by an extended Knowles-Handy lexical scheme.
**
*/
void reorder_ras2(int *docc_in, int *socc_in, int *frozen_docc_in, 
      int *frozen_uocc_in, int *order, int *orbs_per_irrep, 
      int *ras1, int *ras2, int *ras3, int *ras4, int parsed_ras1,
      int parsed_ras2, int do_ras4, int nirreps)
{

   int cnt=0, irrep, point, tmpi;
   int *used, *offset; 
   int *docc, *socc, *frozen_docc, *frozen_uocc;
   int *tras1, *tras2, *tras3, *tras4;

   used = init_int_array(nirreps);
   offset = init_int_array(nirreps);

   docc = init_int_array(nirreps);
   socc = init_int_array(nirreps);
   frozen_docc = init_int_array(nirreps);
   frozen_uocc = init_int_array(nirreps);
   tras1 = init_int_array(nirreps);
   tras2 = init_int_array(nirreps);
   tras3 = init_int_array(nirreps);
   tras4 = init_int_array(nirreps);

   for (irrep=0; irrep<nirreps; irrep++) {
      docc[irrep] = docc_in[irrep];
      socc[irrep] = socc_in[irrep];
      frozen_docc[irrep] = frozen_docc_in[irrep];
      frozen_uocc[irrep] = frozen_uocc_in[irrep];
      }

   /* if the user has not specified RAS I, we must deduce it.
      RAS I does not include any frozen orbitals. */
   if (!parsed_ras1) {
     for (irrep=0; irrep<nirreps; irrep++) {
       ras1[irrep] = docc[irrep] + socc[irrep];
       ras1[irrep] -= frozen_docc[irrep];
     }
   }

   /* if RAS II isn't specified, assume there's not one */
   if (!parsed_ras2) {
     for (irrep=0; irrep<nirreps; irrep++) {
       ras2[irrep] = 0;
     }
   }
   
   /* set up the RAS III or IV array: if RAS IV is used, RAS III must
      be specified and then RAS IV is deduced. */
   for (irrep=0; irrep<nirreps; irrep++) { 
      tmpi = frozen_docc[irrep] + frozen_uocc[irrep] + 
             ras1[irrep] + ras2[irrep];
      if (do_ras4) tmpi += ras3[irrep];
      if (tmpi > orbs_per_irrep[irrep]) {
         fprintf(stderr, "(reorder_ras): orbitals don't add up for irrep %d\n",
            irrep);
         return;
         }
      if (do_ras4) ras4[irrep] = orbs_per_irrep[irrep] - tmpi;
      else ras3[irrep] = orbs_per_irrep[irrep] - tmpi;
      } 

   /* copy ras arrays to tras */
   for (irrep=0; irrep<nirreps; irrep++) {
     tras1[irrep] = ras1[irrep];
     tras2[irrep] = ras2[irrep];
     tras3[irrep] = ras3[irrep];
     tras4[irrep] = ras4[irrep];
   }
   
   /* construct the offset array */
   offset[0] = 0;
   for (irrep=1; irrep<nirreps; irrep++) {
      offset[irrep] = offset[irrep-1] + orbs_per_irrep[irrep-1];
      }
   
   /* do the frozen core */
   for (irrep=0; irrep<nirreps; irrep++) { 
      while (frozen_docc[irrep]) {
         point = used[irrep] + offset[irrep];
         order[point] = cnt++;
         used[irrep]++;
         frozen_docc[irrep]--;
         docc[irrep]--;
         }
      }

   /* do RAS 1: go irrep by irrep */
   for (irrep=0; irrep<nirreps; irrep++) {
     while (tras1[irrep]) {
       point = used[irrep] + offset[irrep];
       order[point] = cnt++;
       used[irrep]++;
       tras1[irrep]--;
     }
   }

   /* do RAS 2: irrep by irrep */
   for (irrep=0; irrep<nirreps; irrep++) {
     while (tras2[irrep]) {
       point = used[irrep] + offset[irrep];
       order[point] = cnt++;
       used[irrep]++;
       tras2[irrep]--;
     }
   }

   /* do RAS 3: irrep by irrep */
   for (irrep=0; irrep<nirreps; irrep++) {
     while (tras3[irrep]) {
       point = used[irrep] + offset[irrep];
       order[point] = cnt++;
       used[irrep]++;
       tras3[irrep]--;
     }
   }

   if (do_ras4) {
     /* do RAS 4: irrep by irrep */
     for (irrep=0; irrep<nirreps; irrep++) {
       while (tras4[irrep]) {
         point = used[irrep] + offset[irrep];
         order[point] = cnt++;
         used[irrep]++;
         tras4[irrep]--;
       }
     }
   }

   /* do frozen uocc */
   for (irrep=0; irrep<nirreps; irrep++) {
     while (frozen_uocc[irrep]) {
       point = used[irrep] + offset[irrep];
       order[point] = cnt++;
       used[irrep]++;
       frozen_uocc[irrep]--;
     }
   }


   /* do a final check */
   for (irrep=0; irrep<nirreps; irrep++) {
     if (used[irrep] > orbs_per_irrep[irrep]) {
       fprintf(stderr, "(reorder_ras): on final check, used more orbitals");
       fprintf(stderr, "   than were available (%d vs %d) for irrep %d\n",
               used[irrep], orbs_per_irrep[irrep], irrep);
     } 
     if (used[irrep] < orbs_per_irrep[irrep]) {
       fprintf(stderr, "(reorder_ras): on final check, used fewer orbitals");
       fprintf(stderr, "   than were available (%d vs %d) for irrep %d\n",
               used[irrep], orbs_per_irrep[irrep], irrep);
     }
   }

   free(used);  free(offset);
   free(docc);  free(socc);  free(frozen_docc);  free(frozen_uocc);
   free(tras1); free(tras2); free(tras3); free(tras4);
}


