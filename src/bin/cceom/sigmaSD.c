
void FSD(int i, int irrep);
void WamefSD(int i, int irrep);
void WmnieSD(int i, int irrep);

/* This function computes the H-bar singles-doubles block contribution
to a Sigma vector stored at Sigma plus 'i' */

void sigmaSD(int i, int irrep) {

  FSD(i, irrep);
  WamefSD(i, irrep);
  WmnieSD(i, irrep);

  return;
}
