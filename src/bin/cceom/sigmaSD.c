
void FSD(int i, int C_irr);
void WamefSD(int i, int C_irr);
void WmnieSD(int i, int C_irr);

/* This function computes the H-bar singles-doubles block contribution
to a Sigma vector stored at Sigma plus 'i' */

void sigmaSD(int i, int C_irr) {

  FSD(i, C_irr);
  WamefSD(i, C_irr);
  WmnieSD(i, C_irr);

  return;
}
