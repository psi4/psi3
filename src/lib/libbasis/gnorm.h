
#ifndef _psi_src_lib_libbasis_gnorm_h_
#define _psi_src_lib_libbasis_gnorm_h_

#include <psitypes.h>

class GaussianNormalization {

  int maxam_;
  PSI_FLOAT* df_;
  PSI_FLOAT** norm_coeffs_;

  // No default constructor
  GaussianNormalization();
  // No assignment operator
  GaussianNormalization& operator=(const GaussianNormalization&);

 public:
  GaussianNormalization(int am);
  ~GaussianNormalization();

  PSI_FLOAT norm(int x, int y, int z);
  PSI_FLOAT norm(int l, int bf);

};

#endif
