
void FDD(int i, int irrep);
void WabefDD(int i, int irrep);
void WmnijDD(int i, int irrep);
void WmbejDD(int i, int irrep);
void WmnefDD(int i, int irrep);

/* This function computes the H-bar doubles-doubles block contribution
to a Sigma vector stored at Sigma plus 'i' */

void sigmaDD(int i, int irrep) {

  FDD(i, irrep);
  WmnijDD(i, irrep);
  WabefDD(i, irrep);
  WmbejDD(i, irrep);
  WmnefDD(i, irrep);

  return;
}
