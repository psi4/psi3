double LHX2Y2(char *cart, int irrep, double omega)
{
  dpdbuf4 X2, Y2, I, W, Z, Z1, Z2, W1, W2;

  sprintf(lbl, "X_%1s_IbjA (-%5.3f)", lbl);
  dpd_buf4_init(&Y2, CC_LR, irrep, 10, 10, 10, 10, 0, lbl);

  dpd_buf4_close(&Y2);

}
