/*
** PARAMS.H
** 
** C. David Sherrill
** Georgia Institute of Technology
** 2001
*/


#define PARM_OUTFILE_MAX           132

/*
** parameters structure: holds user-specified parameters
*/
struct Params {
   char ofname[PARM_OUTFILE_MAX+1]; /* output file name                     */
   char *wfn;               /* wavefunction name                            */
   int print_lvl;           /* print verbosity level                        */ 
   int print_mos;           /* print the molecular orbitals ?               */
   int h_fzc_file;          /* filenum for frozen core operator             */
   int oei_erase;           /* if 1, erase the frozen core operator file    */
   int fzc;                 /* do an implicit frozen core                   */
   int ras_type;            /* define ras I to include 0 or excluded 1 socc */
   /* the next two entries allow for mixing the frozen core fock matrix 
      with the regular fock matrix.  I assume that the regular fock matrix
      is diagonal and given by the orbital eigenvalues; this makes it a 
      lot easier for me to test right now!  I envision setting 
      fzc_fock_coeff and fock_coeff to numbers between 0 and 1 and having
      their sum equal 1.
   */
   double fzc_fock_coeff;   /* coefficient for the frozen core fock matrix  */
   double fock_coeff;       /* coefficient of the regular fock matrix       */
   int mp2nos;              /* if true, get ump2nos                         */
   int unos;                /* if true, get unos                            */
   int canonical;           /* if true, rediag Fock mat in each RAS subspc  */
   int ivo;                 /* get improved virtual orbitals not MVO's      */
   int *docc_virt;          /* treat this doubly occupied orbital as a part */
                            /* of the virtual space during mvo calculations */
  };

