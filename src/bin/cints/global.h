/*-----------------
  global variables
 -----------------*/

#include "data_structs.h"

#ifdef EXTERN
# undef EXTERN
# define EXTERN extern
#else
# define EXTERN
#endif

/*--- Expensive/frequently-used quantities ---*/
EXTERN int ioff[IOFFMAX];             /* Offset array */
EXTERN double fac[CINTS_MAX_AM*2];          /* Factorials */
EXTERN int bc[CINTS_MAX_AM+1][CINTS_MAX_AM+1];    /* Binomial coefficients */
EXTERN double df[MAXFACT*2];          /* Double factorials */
EXTERN int num_ser[CINTS_MAX_AM+2];

/*--- user-specified options ---*/
EXTERN UserOptions_t UserOptions;

/*--- I/O descriptors ---*/
EXTERN IOUnits_t IOUnits;

/*--- Input and Output files ---*/
EXTERN FILE *infile;
EXTERN FILE *outfile;

/*--- Molecule Info ---*/
EXTERN Molecule_t Molecule;

/*--- Basis Set Info ---*/
EXTERN BasisSet_t BasisSet;

/*--- Symmetry arrays ---*/
EXTERN SymmetryInfo_t Symmetry;

/*--- GTO "constant" data (normalization factors, cart2puresphharm coefficients) ---*/
EXTERN GTOs_t GTOs;

/*--- "Taylor"-type evaluator of incomplete gamma function ---*/
EXTERN Fm_Eval_t Taylor_Fm_Eval;

/*--- Information about MOs ---*/
EXTERN MOInfo_t MOInfo;

/*--- DFT information ---*/
EXTERN DFT_options_t DFT_options;

/*--- Calculation-specific 2-index quantities ---*/
EXTERN double **Dens;
EXTERN double **Denso;
EXTERN double **Cocc;     /* Occupied Eigenvector Matrix in AO */
EXTERN double **Cocco;
EXTERN double **Cocca;
EXTERN double **Coccb;
EXTERN double **Lagr;    /* Energy-weighted density of lagrangian in AO basis */
EXTERN double ***ShDens; /* MO shell density */
EXTERN double **G;       /* G-matrix (t.e. part of the Fock matrix) in AO basis */
EXTERN double **Go;      /* open-shell G-matrix in AO basis */
EXTERN double **Grad;    /* Nuclear forces */



