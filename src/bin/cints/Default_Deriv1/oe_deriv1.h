
void oe_deriv(FILE *fpo, double **grad, double **Dens, double **WDens, int nirreps,
	      int num_shells, struct shell_def *shells, struct shell_pair **shell_pairs, struct gaussian_function *cgtos,
	      int num_atoms, struct coordinates *centers, double **bf_norm, int max_am);
double s_ovlp(double a1, double a2, double norm1, double norm2, double ab2);
double overlap_int(double a1, int l1, int m1, int n1, double norm1,
		   double a2, int l2, int m2, int n2, double norm2,
		   struct coordinates AB,
		   struct coordinates PA,
		   struct coordinates PB);
double ke_int(double a1, int l1, int m1, int n1, double norm1,
	      double a2, int l2, int m2, int n2, double norm2,
	      struct coordinates AB,
	      struct coordinates PA,
	      struct coordinates PB);
double f_n(int k, int l1, int l2, double A, double B);
double int_pow(double a, int p);
double ***init_box(int a, int b, int c);
void free_box(double ***box, int a, int b);
