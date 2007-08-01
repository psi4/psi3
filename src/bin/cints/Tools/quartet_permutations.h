/*! \file quartet_permutations.h
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/

#ifdef __cplusplus
extern "C" {
#endif

/** Permute bra and ket in s and put into t */
void ijkl_to_klij(double *s, double *t, int nbra, int nket);

/** Permute bra-indices in s and put into t */
void ijkl_to_jikl(double *s, double *t, int ni, int nj, int nket);

/** Permute ket-indices in s and put into t */
void ijkl_to_ijlk(double *s, double *t, int nbra, int nk, int nl);

#ifdef __cplusplus
}
#endif
