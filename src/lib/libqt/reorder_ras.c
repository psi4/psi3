#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>

/*
** reorder_ras()
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
**
*/
void reorder_ras(int *docc_in, int *socc_in, int *frozen_docc_in, 
      int *frozen_uocc_in, int *order, int *orbs_per_irrep, 
      int *ras1, int *ras2, int *ras3, int *ras4, int do_ras4, int nirreps)
{

   int cnt=0, irrep, point, tmpi;
   int *used, *offset; 
   int *docc, *socc, *frozen_docc, *frozen_uocc;

   used = init_int_array(nirreps);
   offset = init_int_array(nirreps);

   docc = init_int_array(nirreps);
   socc = init_int_array(nirreps);
   frozen_docc = init_int_array(nirreps);
   frozen_uocc = init_int_array(nirreps);

   for (irrep=0; irrep<nirreps; irrep++) {
      docc[irrep] = docc_in[irrep];
      socc[irrep] = socc_in[irrep];
      frozen_docc[irrep] = frozen_docc_in[irrep];
      frozen_uocc[irrep] = frozen_uocc_in[irrep];
      }

   /* set up the RAS III or IV array */
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

   /* do RAS 1: first docc, then socc, then uocc */
   for (irrep=0; irrep<nirreps; irrep++) {
      while (docc[irrep] && ras1[irrep]) {
         point = used[irrep] + offset[irrep];
         order[point] = cnt++;
         used[irrep]++;
         docc[irrep]--;
         ras1[irrep]--;
         }
      }
   for (irrep=0; irrep<nirreps; irrep++) {
      while (socc[irrep] && ras1[irrep]) {
         point = used[irrep] + offset[irrep];
         order[point] = cnt++;
         used[irrep]++;
         socc[irrep]--;
         ras1[irrep]--;
         }
      }
   for (irrep=0; irrep<nirreps; irrep++) {
      while (ras1[irrep]) {
         point = used[irrep] + offset[irrep];
         order[point] = cnt++;
         used[irrep]++;
         ras1[irrep]--;
         }
      }

   /* do RAS 2: first docc, then socc, then uocc */
   for (irrep=0; irrep<nirreps; irrep++) {
      while (docc[irrep] && ras2[irrep]) {
         point = used[irrep] + offset[irrep];
         order[point] = cnt++;
         used[irrep]++;
         docc[irrep]--;
         ras2[irrep]--;
         }
      }
   for (irrep=0; irrep<nirreps; irrep++) {
      while (socc[irrep] && ras2[irrep]) {
         point = used[irrep] + offset[irrep];
         order[point] = cnt++;
         used[irrep]++;
         socc[irrep]--;
         ras2[irrep]--;
         }
      }
   for (irrep=0; irrep<nirreps; irrep++) {
      while (ras2[irrep]) {
         point = used[irrep] + offset[irrep];
         order[point] = cnt++;
         used[irrep]++;
         ras2[irrep]--;
         }
      }

   /* do RAS 3: first docc, then socc, then uocc */
   for (irrep=0; irrep<nirreps; irrep++) {
      while (docc[irrep] && ras3[irrep]) {
         point = used[irrep] + offset[irrep];
         order[point] = cnt++;
         used[irrep]++;
         docc[irrep]--;
         ras3[irrep]--;
         }
      }
   for (irrep=0; irrep<nirreps; irrep++) {
      while (socc[irrep] && ras3[irrep]) {
         point = used[irrep] + offset[irrep];
         order[point] = cnt++;
         used[irrep]++;
         socc[irrep]--;
         ras3[irrep]--;
         }
      }
   for (irrep=0; irrep<nirreps; irrep++) {
      while (ras3[irrep]) {
         point = used[irrep] + offset[irrep];
         order[point] = cnt++;
         used[irrep]++;
         ras3[irrep]--;
         }
      }

   if (do_ras4) {
      /* do RAS 4: first docc, then socc, then uocc */
      for (irrep=0; irrep<nirreps; irrep++) {
         while (docc[irrep] && ras4[irrep]) {
            point = used[irrep] + offset[irrep];
            order[point] = cnt++;
            used[irrep]++;
            docc[irrep]--;
            ras4[irrep]--;
            }
         }
      for (irrep=0; irrep<nirreps; irrep++) {
         while (socc[irrep] && ras4[irrep]) {
            point = used[irrep] + offset[irrep];
            order[point] = cnt++;
            used[irrep]++;
            socc[irrep]--;
            ras4[irrep]--;
            }
         }
      for (irrep=0; irrep<nirreps; irrep++) {
         while (ras4[irrep]) {
            point = used[irrep] + offset[irrep];
            order[point] = cnt++;
            used[irrep]++;
            ras4[irrep]--;
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
}


