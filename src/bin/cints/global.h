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
EXTERN double fac[MAX_AM*2];          /* Factorials */
EXTERN int bc[MAX_AM+1][MAX_AM+1];    /* Binomial coefficients */
EXTERN double df[MAXFACT*2];          /* Double factorials */
EXTERN int num_ser[MAX_AM+2];

/*--- user-specified options ---*/
EXTERN UserOptions_t UserOptions;

/*--- I/O descriptors ---*/
EXTERN IOUnits_t IOUnits;

/*--- Input and Output files ---*/
EXTERN FILE *infile;
EXTERN FILE *outfile;

/*--- used to pass data to functions in LIBINT ---*/
EXTERN double *int_stack;
EXTERN double *zero_stack;             /* LIBDERIV only */
EXTERN double AB[3], CD[3];
EXTERN prim_data *Shell_Data;
EXTERN double *ABCD[12];               /* LIBDERIV only */
EXTERN double *te_ptr[NUM_TE_TYPES];   /* LIBR12 only */
EXTERN double AB[3], CD[3], AC[3];     /* LIBR12 only */
EXTERN double ABdotAC, CDdotCA;        /* LIBR12 only */

/*--- Table of Fj(T) values ---*/
EXTERN double_array_t int_fjttable;

/*--- Molecule Info ---*/
EXTERN Molecule_t Molecule;

/*--- Basis Set Info ---*/
EXTERN BasisSet_t BasisSet;

/*--- Symmetry arrays ---*/
EXTERN SymmetryInfo_t Symmetry;

/*--- GTO "constant" data (normalization factors, cart2puresphharm coefficients) ---*/
EXTERN GTOs_t GTOs;

/*--- Information about MOs ---*/
EXTERN MOInfo_t MOInfo;

/*--- Calculation-specific 2-index quantities ---*/
EXTERN double **Dens;    /* Total density in AO basis */
EXTERN double **Denso;   /* Open-shell density in AO basis */
EXTERN double **Lagr;    /* Energy-weighted density of lagrangian in AO basis */
EXTERN double ***ShDens; /* MO shell density */
EXTERN double **G;       /* in SO basis */
EXTERN double **Go;      /* in SO basis */
EXTERN double **Grad;    /* Nuclear forces */
