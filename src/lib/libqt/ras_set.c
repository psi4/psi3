/*!
  \file ras_set.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>

/*!
** ras_set()
**
** This function sets up the number of orbitals per irrep for each of the
** RAS subspaces [frozen core, RAS I, RAS II, RAS III, RAS IV, frozen virts].
** It also obtains the appropriate orbital reordering array.  The 
** reordering array takes a basis function in Pitzer ordering (orbitals 
** grouped according to irrep) and gives the corresponding index 
** in the RAS numbering scheme.  Orbitals are numbered according to 
** irrep within each of the subspaces.
**
** C. David Sherrill
** Center for Computational Quantum Chemistry
** University of Georgia, 25 June 1995
**
** Parameters:
**  \param nirreps     =  num of irreps in computational point group
**  \param nbfso       =  num of basis functions in symmetry orbitals (num MOs) 
**  \param freeze_core =  1 to remove frozen core orbitals from ras_opi
**  \param orbspi      =  array giving num symmetry orbitals (or MOs) per irrep
**  \param docc        =  array of doubly occupied orbitals per irrep
**  \param socc        =  array of singly occupied orbitals per irrep
**  \param frdocc      =  array of frozen core per irrep
**  \param fruocc      =  array of frozen virtuals per irrep
**  \param ras_opi     =  matrix giving num of orbitals per irrep per ras space,
**                        addressed as ras_opi[ras_space][irrep]
**  \param order       =  array nbfso big which maps Pitzer to Correlated order
**  \param ras_type    =  if 1, put docc and socc together in same RAS space 
**                        (RAS I), as appropriate for DETCI.  If 0, put socc
**                        in its own RAS space (RAS II), as appropriate for CC.
**
** Returns: 1 for success, 0 otherwise
**
*/
int ras_set(int nirreps, int nbfso, int freeze_core, int *orbspi,
      int *docc, int *socc, int *frdocc, int *fruocc, 
      int **ras_opi, int *order, int ras_type)
{
  int i, irrep, point, tmpi, cnt=0;
  int errcod, errbad, parsed_ras1=0, parsed_ras2=0, do_ras4;
  int *used, *offset, **tras;

  used = init_int_array(nirreps);
  offset = init_int_array(nirreps);
  
  /* zero everything out, don't assume the user has done it already */
  zero_int_array(docc, nirreps);
  zero_int_array(socc, nirreps);
  zero_int_array(frdocc, nirreps);
  zero_int_array(fruocc, nirreps);
  for (i=0; i<4; i++) {
    zero_int_array(ras_opi[i], nirreps);
  }
  zero_int_array(order, nbfso);
  
  
  /* now use the parser to get the arrays we require */
  errcod = ip_int_array("DOCC",docc,nirreps); 
  if (errcod != IPE_OK) {
    fprintf(stderr, "(ras_set): Error reading DOCC\n");
    exit(1);
  }
  errcod = ip_int_array("SOCC",socc,nirreps);
  errcod = ip_int_array("FROZEN_DOCC",frdocc,nirreps);
  errcod = ip_int_array("FROZEN_UOCC",fruocc,nirreps);
  
  errbad=0; do_ras4=1;
  errcod = ip_int_array("RAS1", ras_opi[0], nirreps);
  if (errcod == IPE_OK) parsed_ras1 = 1;
  else if (errcod == IPE_KEY_NOT_FOUND) {
    errcod = ip_int_array("ACTIVE_CORE", ras_opi[0], nirreps);
    if (errcod == IPE_OK) parsed_ras1 = 1;
    else if (errcod != IPE_KEY_NOT_FOUND) errbad = 1;
  }
  else errbad = 1;
  errcod = ip_int_array("RAS2", ras_opi[1], nirreps);
  if (errcod == IPE_OK) parsed_ras2 = 1;
  else if (errcod == IPE_KEY_NOT_FOUND) {
    errcod = ip_int_array("MODEL_SPACE", ras_opi[1], nirreps);
    if (errcod == IPE_OK) parsed_ras2 = 1;
    else if (errcod != IPE_KEY_NOT_FOUND) errbad = 1;
  }
  else errbad = 1;
  errcod = ip_int_array("RAS3", ras_opi[2], nirreps);
  if (errcod != IPE_OK && errcod != IPE_KEY_NOT_FOUND) errbad=1;
  if (errcod == IPE_KEY_NOT_FOUND) do_ras4=0;
  
  if (errbad == 1) {
    fprintf(stderr, "(ras_set): trouble parsing RAS keyword\n");
    return(0);
  }
  
  /* if the user has not specified RAS I, we must deduce it.
   * RAS I does not include any FZC orbs but does include COR orbs
   */
  
  if (!parsed_ras1) {
    for (irrep=0; irrep<nirreps; irrep++) {
      if (ras_type==1) ras_opi[0][irrep] = docc[irrep] + socc[irrep];
      if (ras_type==0) ras_opi[0][irrep] = docc[irrep];
      ras_opi[0][irrep] -= frdocc[irrep]; /* add back later for COR */ 
    }
  }
  
  
  /* if the user hasn't specified RAS II, look for val_orb             */
  /* val_orb should typically be RAS I + RAS II, so subtract out RAS I */ 
  if (!parsed_ras2) {
    errcod = ip_int_array("VAL_ORB",ras_opi[1],nirreps);
    if (errcod != IPE_OK) {
      for (irrep=0; irrep<nirreps; irrep++) {
	if (ras_type==1) ras_opi[1][irrep] = 0;
	if (ras_type==0) ras_opi[1][irrep] = socc[irrep];
      }
    }
    else {
      for (irrep = 0; irrep<nirreps; irrep++) { 
         ras_opi[1][irrep] -= ras_opi[0][irrep];
         if (ras_opi[1][irrep] < 0) {
           fprintf(stderr, "(ras_set): val_orb must be larger than RAS I\n");
           return(0);
         }
      }
    }
  }
  
  /* set up the RAS III or IV array: if RAS IV is used, RAS III must
   * be specified and then RAS IV is deduced.
   */
  
  for (irrep=0; irrep<nirreps; irrep++) { 
    tmpi = frdocc[irrep] + fruocc[irrep] + ras_opi[0][irrep] + 
      ras_opi[1][irrep];
    if (do_ras4) tmpi += ras_opi[2][irrep];
    if (tmpi > orbspi[irrep]) {
      fprintf(stderr, "(ras_set): orbitals don't add up for irrep %d\n",
	      irrep);
      return(0);
    }
    if (do_ras4) ras_opi[3][irrep] = orbspi[irrep] - tmpi;
    else ras_opi[2][irrep] = orbspi[irrep] - tmpi;
  } 
  
  /* copy everything to the temporary RAS arrays: */
  /* add subspaces for frozen orbitals            */
  tras = init_int_matrix(6, nirreps);
  for (irrep=0; irrep<nirreps; irrep++) {
    tras[0][irrep] = frdocc[irrep];
    tras[5][irrep] = fruocc[irrep];
  }
  for (i=0; i<4; i++) {
    for (irrep=0; irrep<nirreps; irrep++) {
      tras[i+1][irrep] = ras_opi[i][irrep];
    }
  }
 
  /* construct the offset array */
  offset[0] = 0;
  for (irrep=1; irrep<nirreps; irrep++) {
    offset[irrep] = offset[irrep-1] + orbspi[irrep-1];
  }
  
  for (i=0; i<6; i++) {
    for (irrep=0; irrep<nirreps; irrep++) { 
      while (tras[i][irrep]) {
	point = used[irrep] + offset[irrep];
	if (point < 0 || point >= nbfso) {
	  fprintf(stderr, "(ras_set): Invalid point value\n");
	  exit(1);
	}
	order[point] = cnt++;
	used[irrep]++;
	tras[i][irrep]--;
      }
    }
  }
  

  /* do a final check */

  for (irrep=0; irrep<nirreps; irrep++) {
    if (used[irrep] > orbspi[irrep]) {
      fprintf(stderr, "(ras_set): on final check, used more orbitals");
      fprintf(stderr, "   than were available (%d vs %d) for irrep %d\n",
	      used[irrep], orbspi[irrep], irrep);
      errbad = 1;
    } 
    if (used[irrep] < orbspi[irrep]) {
      fprintf(stderr, "(ras_set): on final check, used fewer orbitals");
      fprintf(stderr, "   than were available (%d vs %d) for irrep %d\n",
	      used[irrep], orbspi[irrep], irrep);
      errbad = 1;
    }
  }
  
  /* for restricted COR orbitals */
  if (!freeze_core) {
    for (irrep=0; irrep<nirreps; irrep++) {
      ras_opi[0][irrep] += frdocc[irrep];
    }
  }
  
  free(used);  free(offset);
  free_int_matrix(tras, 6);

  return(!errbad);

}


