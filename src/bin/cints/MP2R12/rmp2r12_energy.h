/*! \file rmp2r12_energy.h
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#define LOCK_RS_SHELL 0      /*--- When updating (js|ia) and (jr|ia) in (JS|IA) lock blocks corresponding to the entire shell
			       blocks or just the appropriate basis functions (more fine-grained in the second case) ---*/

void rmp2r12_energy();

typedef struct {
    int num_i_per_ibatch;
    int num_ibatch;
    int ibatch_first;
    int num_arrived;
} RMP2R12_Status_t;
