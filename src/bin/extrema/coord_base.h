/*#############################################################################
  coord_base.h

  data and functions common to all coordinate systems
  ###########################################################################*/

double **symm_matrix_invert(double**,int,int,int);

class coord_base {
  
  protected:
    
    double *coords, *carts, *grads, *c_grads, **H, **H_evec, *H_eval, **u, *coords_old, *grads_old, **H_old, *coord_write, *coord_temp, *masses;

  public:

    coord_base();
    ~coord_base();    
    void print_carts(double conv);
    void print_c_grads();
    void print_grads();
    void diagonalize_H();
    void print_H();
    void print_u();
    void read_carts();
    void read_opt();
    void write_opt();
    virtual void opt_step();
    virtual double* compute_s();
    void update_bfgs();
    void update_ms();
    virtual void read_file11();
    virtual void write_file30();
    void grad_test();
};












