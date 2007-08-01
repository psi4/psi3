/*! \file Params.h
    \ingroup (MP2R12)
    \brief Enter brief description of file here 
*/
/* Struct for input parameters */
struct Params {
    int print_lvl;          /* Printing level */
    int c_limit;            /* Whether to use the limiting form for the B^{-1}*V */
    double tolerance;       /* Cutoff for reading in integrals */
    char *wfn;              /* The wavefunction */
    int keep_integrals;     /* Keep the integrals? */
};
