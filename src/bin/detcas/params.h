/*
** PARAMS.H
** 
** C. David Sherrill
** University of California, Berkeley
** 1998
*/


#define PARM_OUTFILE_MAX           132

/*
** parameters structure: holds user-specified parameters
*/
struct params {
   char ofname[PARM_OUTFILE_MAX+1]; /* output file name                     */
   char *dertype;           /* derivative level: none, first, etc.          */
   int print_lvl;           /* print verbosity level                        */ 
   int print_mos;           /* print the molecular orbitals ?               */
   int rms_grad_convergence;/* convergence, 10^-n, on RMS of orbital grad   */
   int energy_convergence;  /* convergence, 10^-n, on CI energy             */
   int oei_file;            /* file number for one-electron integrals       */
   int oei_erase;           /* erase onel ints after reading them?          */
   int tei_file;            /* file number for two-electron integrals       */
   int tei_erase;           /* erase twoel ints after reading them?         */
   int opdm_file;           /* file number for one-particle density matrix  */
   int opdm_erase;          /* erase onepdm ints after reading?             */
   int tpdm_file;           /* file number for two-particle density matrix  */
   int tpdm_erase;          /* erase twopdm after reading?                  */
   int lag_file;            /* file number for lagrangian                   */
   int lag_erase;           /* erase lagrangian after reading?              */
   int fci;                 /* do a FULL ci calc?  (affects independent prs */
   int fzc;                 /* do implicit frozen core (remove those orbs)? */
                            /* the alternative is a "restricted core" calc  */
   int filter_ints;         /* filter out the frozen orbital integrals?     */
   int scale_grad;          /* scale the orbital gradient by the appx Hess? */
   int diis_start;          /* how many diis vectors built up before start  */
   int diis_freq;           /* how many iters to go before a diis step      */
   int diis_min_vecs;       /* how many vectors required before do diis?    */
   int diis_max_vecs;       /* how many vectors maximum to hold?            */
   double scale_step;       /* stepsize scaling factor                      */
  };

