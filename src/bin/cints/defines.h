/*----------------------------------
  Flags for conditional compilation
 ----------------------------------*/
#define PRINT 0          /* Print two-electron integrals out? */
#define USE_MM 1         /* Use matrix multiplies rather than just nested loops
                            in cart->puream transformation and normalization (should always be
			    set to 1 except for testing purposes)*/
#define SPARSE_C2P 1     /* Use sparcity of cartesian to spherical harmonics transformation */
#define USE_BLAS 1       /* Use routines from vendor-provided BLAS library for
			    4-index transformations such as cart->puream
			    (only if SPARSE_C2P is set 0, otherwise sparse matrix
			    multiplies will be used), AO->MO, etc.;
			    Use only if you have libblas.a or its other analog available */
#define SCF_ONLY 1       /* If you want to be able to compute integrals needed in
			    SCF only if WFN=SCF - set this to 1 */

/*---------------------
  Predefined constants
 ---------------------*/
#define ZERO 1E-15                 /* Definition of a floating-point "zero" */
#define EPS 1.0e-17                /* Another definition of floating-point "zero"
				      used in computing auxiliary function */
#define CUTOFF 15                  /* Default cutoff on the integrals */
#define ROT_INV_TOLER 1E-4         /* Tolerance on the "rotational variance" */

/*----------------------------------
  Thresholds for printing out stuff
 ----------------------------------*/
#define PRINT_INTRO 1    /* Print level to print intro overhead */
#define PRINT_OPTIONS 1  /* Print level to print options out */
#define PRINT_BASIS 3    /* Print level to print basis set information */
#define PRINT_GEOMETRY 2 /* Print level to print cartesian geometry */
#define PRINT_CCOEFF 4   /* Print level to print coupling coefficients */
#define PRINT_OPDM 4     /* Print level to print onepdms */
#define PRINT_OEI 3      /* Print level to print one-electron integrals */
#define PRINT_OEDERIV 2      /* Print level to print oe and nuclear contrinution to gradient */
#define PRINT_TEDERIV 2      /* Print level to print te contribution to gradient  */
#define PRINT_MOINFO_CORR 1  /* Print level to print orbital information for correlated calculations */

/*-------------------------
  Default sizes for arrays
 -------------------------*/
#define MAXFACT 100
#define MAXNIRREPS 8
#define MAX_NUM_DOUBLES 2500000 /* Default number of double words to use in MP2 */
#define IOFFMAX 16384

/*----------------
  Macro functions
 ----------------*/
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#define INDEX(a,b) (((a)>(b)) ? ioff[(a)] + (b) : ioff[(b)] + (a))
#define SWAP(a,b) {dum = (a); (a) = (b); (b) = dum;}

/*------
  Misc.
 ------*/
#define NUM_TE_TYPES 4   /* Number of types of integrals to compute
			    0 - ERIs
			    1 - r12
			    2 - [r12,T1]
			    3 - [r12,T2]
			    */
