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
  };

