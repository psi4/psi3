#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_scf_ptrs(): Read the pointers to the data stored in the scf section of file30.
** 
** Here is a map of this section*/

/*  The structure of the SCF section is as follows:
   (listed first is the length in integer words, 
   variable names correspond to those given in CSCF3.0)


   | 20 (label)  | 40 (don't know)  |  20 (pointers written    |
   |             |                  |      this function)      |           
   _____________________________________________________________


   |  6*natoms      | 8 ( 4 double       | 2 (7 char string   |  mxcoef*2           |
   |  (geometry)    | constants Enuc,    | + '\0' to identify |                     |
   |                | Escf, Eref, Ecorr) | corr. energy Ecorr | (Alpha SCF eigenvec)|  
   ----------------------------------------------------------------------------------

   | mxcoef*2           | nmo*2               | nmo*2              |
   | (Beta SCF eigenvec)| (Alpha SCF eigenval)| (Beta SCF eigenval)|
   -----------------------------------------------------------------

   | num_ir*4               | n_so_typs                  |
   | (Non-zero irrep labels)| (number of orbitals/irrep) |
   -------------------------------------------------------
   
   | n_so_typs                   | n_so_typs                 |
   | (number closed shells/irrep)| (number open shells/irrep)|
   -----------------------------------------------------------

   | ioff[n_open]*2 | ioff[n_open]*2 | nmo*(nmo+1)/2*2   |
   | (Alpha coefs)  | (Beta coefs)   | (Alpha Lagrangian)|
   -------------------------------------------------------

   | nmo*(nmo+1)/2*2  |
   | (Beta Lagrangian)|
   --------------------*/

   /*--------------------------------------------------------------*/
   /* Pointers will be an array that will contain the positions in */
   /*   file30 for the following informaiton                       */
   /* [0] = Alpha SCF eigenvector (the only one stored if not UHF) */
   /* [1] = Beta  SCF eigenvector                                  */
   /* [2] = Alpha SCF eigenvalues (the only one stored if not UHF) */
   /* [3] = Beta  SCF eigenvalues                                  */
   /* [4] = Irrep Labels for irreps with non-zero occupancies      */
   /* [5] = Number of orbitals per irrep for non-zero occupancies  */
   /* [6] = Number of closed shells per irrep with non-zero occ    */
   /* [7] = Number of open shells per irrep with non-zero occ      */
   /* [8] = Alpha Coupling Constants                               */
   /* [9] = Beta  Coupling Constants                               */
   /* [10] = Alpha Lagrangian (the only one stored if not UHF)     */
   /* [11] = Beta  Lagrangian                                      */
   /* It is important to note that the pointers stored in the file */
   /* will not be absolute pointer values in the C sense.  They    */
   /* will actually be the pointer value/sizeof(int) +1 to keep    */
   /* consistent with FORTRAN.  So when the pointers need be used  */
   /* they need to be converted to the correct value by the formula*/
   /*           value = (value from file30)*sizeof(int) -1         */
   /*--------------------------------------------------------------*/
  
/*   returns an array containing the pointers */

int *file30_rd_scf_ptrs(void)
{
    PSI_FPTR scf_point, junk;
    int *scf_ptrs;
    
    scf_ptrs = (int *)malloc(sizeof(int)*20);
    
    scf_point = (PSI_FPTR) (info30_.mcalcs[0]+60-1)*sizeof(int);
    
    wreadw(info30_.filenum,(char *) scf_ptrs, sizeof(int)*20,scf_point, &junk);
    
    return scf_ptrs;
}
    

