#ifndef _psi_src_bin_cints_MP2R12_rmp2r12_energy_h
#define _psi_src_bin_cints_MP2R12_rmp2r12_energy_h

/*! \file rmp2r12_energy.h
    \ingroup (CINTS)
*/
#define LOCK_RS_SHELL 0      /*--- When updating (js|ia) and (jr|ia) in (JS|IA) lock blocks corresponding to the entire shell
			       blocks or just the appropriate basis functions (more fine-grained in the second case) ---*/
namespace psi { namespace CINTS {
void rmp2r12_energy();

typedef struct {
    int num_i_per_ibatch;
    int num_ibatch;
    int ibatch_first;
    int num_arrived;
} RMP2R12_Status_t;
};}
#endif
