
void setup();
void start_io(int argc, char *argv[]);
void stop_io();
void punt(char* message);
double distance_calc(struct coordinates g1, struct coordinates g2);
double ***init_box(int, int, int);
void free_box(double ***, int, int);
