/* $Id$ */
/* $Log$
 * Revision 1.2  2002/05/10 05:44:06  crawdad
 * Changed to variable name "nint" to "nnint" to avoid conflict with typedef
 * in Tru64's math.h.
 * -TDC
 *
/* Revision 1.1.1.1  2000/02/04 22:51:33  evaleev
/* Started PSI 3 repository
/*
/* Revision 1.1  1991/06/15 22:06:42  seidl
/* Initial revision
/* */

#define MAX_BASIS 200
#define MAX_STRING 512

#ifdef EXTERN
# undef EXTERN
# define EXTERN extern
#else
# define EXTERN
#endif

EXTERN FILE *infile, *outfile;

EXTERN double repnuc;           /* nuclear repulsion */

EXTERN char blabel[81];

EXTERN int print;               /* print flag */
EXTERN int toler;               /* max iterations */

EXTERN int n_so_typs;           /* number of irreps w/ non-zero num of so's */
EXTERN int num_ir;              /* # of symmetry types */
EXTERN int pos34;               /* pointer to location in file34 */
EXTERN int nnint;                /* number of pki ints in present batch */
EXTERN int nbfao;
EXTERN int nbfso;
EXTERN int nbatri;
EXTERN int nbstri;
EXTERN int nbasis,ntri;
EXTERN int maxbuf;
EXTERN int mxcoef;
EXTERN int nsect;
EXTERN int natom;
EXTERN int iopen;
EXTERN int twocon;
EXTERN int ntypes;
EXTERN int hsos;
EXTERN int ci_calc;

EXTERN int itap30,itap40,itap34,itap37;              

EXTERN int ioff[1024];          /* matrix offsets */
EXTERN int degen[20];           /* degeneracy of each irrep */
EXTERN int num_so[20];          /* # so's per irrep */
EXTERN int block_num[20];
EXTERN int ideg[10];            /* offsets for each irrep */
EXTERN int nlamda[10],nclosd[10],nopen[10];
EXTERN int nsorb[10];
EXTERN int *block_locs;

EXTERN double enuc,escf;

EXTERN double *e_vals,*occ_num;
EXTERN double **e_vecs_so,**e_vecs_ao,**ao_to_so;
EXTERN double *charges,**coords;
EXTERN double *alpha,*beta,*alpb,*betb;
EXTERN double **zeta_so,**densa,**densb;
EXTERN double two_occ[10];
EXTERN double alpc[15],betc[15];
EXTERN double focc[10];


EXTERN union psi_buffer {
          int *lbli;
          double *pki;
          } oubuf;
