/*****************************************************************************

	CARTESIANS.H

        class declaration for cartesian
******************************************************************************/

class cartesians {
    double energy;
    int num_atoms;
    double *atomic_num;
    double *coord;
    double *grad;
    double *mass;

  public:
    ~cartesians() {free(atomic_num); free(coord); free(grad); free(mass); }
    void print(int flag, FILE *fp_out, int new_geom_file, char *disp_label,
               int disp_num);
    void set_coord(double *geom) {
      int i;
      for (i=0;i<num_atoms*3;++i)
        coord[i] = geom[i];
      return;
    }
    double *get_coord() {
      int i;
      double *copy;
      copy = init_array(num_atoms*3);
      for (i=0;i<num_atoms*3;++i)
        copy[i] = coord[i];
      return copy;
    }
    double *get_mass() {
      int i;
      double *copy;
      copy = init_array(num_atoms*3);
      for (i=0;i<num_atoms*3;++i)
        copy[i] = mass[i];
      return copy;
    }
    void mult(double factor) {
      int i;
      for (i=0;i<num_atoms*3;++i) {
         coord[i] *= factor;
      }
      return;
    }
    double val(int i, int j) { return coord[3*i+j]; }
    double *forces();
    int get_num_atoms() {return num_atoms; }
    void set_num_atoms(int new_num) {num_atoms = new_num;}
    void set_energy(double new_energy) {energy = new_energy;}
    double get_energy();
    double get_atomic_num(int i) { return atomic_num[i]; }
    cartesians();
};

