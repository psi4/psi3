#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>


/*
** reorder_qt()
**
** This function constructs a reordering array according to the
** "Quantum Trio" standard ordering, in which the orbitals are divided
** into the following sets: frozen core, then doubly occupied, then singly
** occupied, then virtuals, then deleted (frozen) virtuals.
** The reordering array takes a basis function in 
** Pitzer ordering (orbitals grouped according to irrep) and gives the
** corresponding index in the Quantum Trio numbering scheme.
**
** Should give the same reordering array as in the old libread30 routines.
**
** C. David Sherrill
** Center for Computational Quantum Chemistry
** University of Georgia, 1995
**
*/
void reorder_qt(int *docc_in, int *socc_in, int *frozen_docc_in, 
      int *frozen_uocc_in, int *order, int *orbs_per_irrep, int nirreps)
{

   int cnt=0, irrep, point, tmpi;
   int *used, *offset; 
   int *docc, *socc, *frozen_docc, *frozen_uocc; 
   int *uocc;

   used = init_int_array(nirreps);
   offset = init_int_array(nirreps);

   docc = init_int_array(nirreps);
   socc = init_int_array(nirreps);
   frozen_docc = init_int_array(nirreps);
   frozen_uocc = init_int_array(nirreps);
   uocc = init_int_array(nirreps);

   for (irrep=0; irrep<nirreps; irrep++) {
      docc[irrep] = docc_in[irrep];
      socc[irrep] = socc_in[irrep];
      frozen_docc[irrep] = frozen_docc_in[irrep];
      frozen_uocc[irrep] = frozen_uocc_in[irrep];
      }
 
   /* construct the offset array */
   offset[0] = 0;
   for (irrep=1; irrep<nirreps; irrep++) {
      offset[irrep] = offset[irrep-1] + orbs_per_irrep[irrep-1];
      }
   
   /* construct the uocc array */
   for (irrep=0; irrep<nirreps; irrep++) {
      tmpi = frozen_uocc[irrep] + docc[irrep] + socc[irrep];
      if (tmpi > orbs_per_irrep[irrep]) {
         fprintf(stderr, "(reorder_qt): orbitals don't add up for irrep %d\n",
            irrep);
         return;
         }
      else
         uocc[irrep] = orbs_per_irrep[irrep] - tmpi;
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

   /* do doubly occupied orbitals */
   for (irrep=0; irrep<nirreps; irrep++) {
      while (docc[irrep]) {
         point = used[irrep] + offset[irrep];
         order[point] = cnt++;
         used[irrep]++;
         docc[irrep]--;
         }
      }

   /* do singly-occupied orbitals */
   for (irrep=0; irrep<nirreps; irrep++) {
      while (socc[irrep]) {
         point = used[irrep] + offset[irrep];
         order[point] = cnt++;
         used[irrep]++;
         socc[irrep]--;
         }
      }

   /* do virtual orbitals */
   for (irrep=0; irrep<nirreps; irrep++) {
      while (uocc[irrep]) {
         point = used[irrep] + offset[irrep];
         order[point] = cnt++;
         used[irrep]++;
         uocc[irrep]--;
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
         fprintf(stderr, "(reorder_qt): on final check, used more orbitals");
         fprintf(stderr, "   than were available (%d vs %d) for irrep %d\n",
            used[irrep], orbs_per_irrep[irrep], irrep);
         } 
      }

   free(used);  free(offset);
   free(docc);  free(socc);  free(frozen_docc);  free(frozen_uocc);
   free(uocc);
}


