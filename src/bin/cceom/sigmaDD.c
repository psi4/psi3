
void FDD(int i, int C_irr);
void WabefDD(int i, int C_irr);
void WmnijDD(int i, int C_irr);
void WmbejDD(int i, int C_irr);
void WmnefDD(int i, int C_irr);

/* This function computes the H-bar doubles-doubles block contribution
to a Sigma vector stored at Sigma plus 'i' */

void sigmaDD(int i, int C_irr) {

  FDD(i, C_irr);
  WmnijDD(i, C_irr);
  WabefDD(i, C_irr);
  WmbejDD(i, C_irr);
  WmnefDD(i, C_irr);

  return;
}
