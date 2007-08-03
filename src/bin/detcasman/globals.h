/*! \file 
    \ingroup (DETCASMAN)
    \brief Enter brief description of file here 
*/
/*
** GLOBALS.H
**
** List of all the global data used by the program
**
** Note that these are given as "extern", so they must be defined
** in the main program!
**
** C. David Sherrill
** University of California, Berkeley
*/

extern FILE *infile, *outfile;
extern char *psi_file_prefix;
extern int converged; 
extern int ncasiter;          /* max cas iterations */
extern char detci_string[80]; /* string containing system call for DETCI  */
extern double ci_conv;        /* desired CI convergence 
                                 (changes dynamically during CAS opt)     */
extern double scale_conv;     /* CI convergence threshold = 
                                 orbital gradient * scale_conv            */
