#ifndef globals_h
#define globals_h

#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

EXTERN double*** Christoffel();
EXTERN double dot_prod(double *, double *);
EXTERN double* vect_prod(double *, double *);

#endif

