void parsing();
void grid_unitvec();
void compute_density();
void read_density();
void get_nmo();
void compute_overlap();
void compute_onecgt();
void read_basset_info();
void read_zvec();
void MI_OSrecurs(double, double, double, double, double, double, double, int, int, int);
void AI_OSrecurs(double, double, double, double, double, double, double, double, double, double, int, int);
void populate();
void compute_mp_ref_xyz();
void move2ref();
void init_xyz();
double ***init_box(int, int, int);
void free_box(double ***, int, int);

