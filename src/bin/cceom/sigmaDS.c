
void WmaijDS(int i, int C_irr);
void WabejDS(int i, int C_irr);
void WbmfeDS(int i, int C_irr);
void WnmjeDS(int i, int C_irr);

/* This function computes the H-bar doubles-singles block contribution
to a Sigma vector stored at Sigma plus 'i' */

void sigmaDS(int i, int C_irr) {

  WmaijDS(i, C_irr);
  WabejDS(i, C_irr);
  WnmjeDS(i, C_irr);
  WbmfeDS(i, C_irr);

  return;
}
