#ifndef _psi_src_bin_cints_mkpt2_ints_h
#define _psi_src_bin_cints_mkpt2_ints_h

#include <vector>
#include <string>

/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here
*/
#define LOCK_RS_SHELL 0      /*--- When updating (js|ia) and (jr|ia) in (JS|IA) lock blocks corresponding to the entire shell
			       blocks or just the appropriate basis functions (more fine-grained in the second case) ---*/
namespace psi {
  namespace CINTS {
    namespace mkpt2 {
      void correlate(char *ptgrp, int irrep, int& nirreps_old, int& nirreps_new,int*& correlation);
      void read_mo_space(int nirreps_ref, int& n, std::vector<int>& mo, std::string labels);
      std::vector<int> read_chkpt_intvec(int n, int* array);
      bool space(char c);
      bool not_space(char c);
      std::vector<std::string> split(const std::string& str);

      void mkpt2_ints(void);
      typedef struct {
        int num_i_per_ibatch;
        int num_ibatch;
        int num_arrived;
      } MkPT2_Status_t;
    }
  }
}
#endif
