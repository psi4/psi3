
void WmaijDS(int i, int irrep);
void WabejDS(int i, int irrep);
void WbmfeDS(int i, int irrep);
void WnmjeDS(int i, int irrep);

/* This function computes the H-bar doubles-singles block contribution
to a Sigma vector stored at Sigma plus 'i' */

void sigmaDS(int i, int irrep) {

  WmaijDS(i, irrep);
  WabejDS(i, irrep);
  WnmjeDS(i, irrep);
  WbmfeDS(i, irrep);

  return;
}
