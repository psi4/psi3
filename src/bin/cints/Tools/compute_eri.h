#ifndef __cplusplus
#error "Tools/compute_eri.h cannot be used in C programs"
#endif

/// Returns the number of integrals computed
extern "C" int compute_eri(double* target, Libint_t* Libint, int& si, int& sj, int& sk, int& sl,
                           int& inc1, int& inc2, int& inc3, int& inc4, const bool do_not_permute);
