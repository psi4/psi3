/*#############################################################################
  internals.h

  derived internals class;
    that which is common to both zmatrix and delocalized coodinate systems
  ###########################################################################*/

class internals : public coord_base {

  protected:
    double **full_geom, **B, **G, **A;

  public:
    internals() : coord_base() {
	full_geom = init_matrix(num_entries,3);
        B = (double**) malloc(num_coords*sizeof(double*));
	G = init_matrix(num_coords,num_coords);
        A = init_matrix(3*num_entries,num_coords);
	return;
    }
    ~internals() {
	free_matrix(full_geom,num_entries);
        free_matrix(B,num_coords);
	free_matrix(G,num_coords);
	free_matrix(A,3*num_entries);
	return;
    }
    virtual void compute_B() = 0;
    virtual void cart_to_internal(double*) = 0;
    virtual void print_internals() = 0;
    void back_transform();
    double *B_row_bond(double*, int, int);
    double *B_row_angle(double*, int, int, int);
    double *B_row_tors(double*, int, int, int, int);
    void print_B();
    void compute_G();
    void print_G();
    void compute_A();
    void print_A();
    void grad_trans();
    void optimize_internals(internals* icrd);
    virtual void initial_H() = 0;
};






