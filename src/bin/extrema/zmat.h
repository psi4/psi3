/*#############################################################################
  zmat.h

  derived z matrix class
  ###########################################################################*/

class zmat : public internals {

    struct z_entry* z_geom;
    char **felement;
    simple *simples;

  public:

    zmat();
    ~zmat() {
        int i; 
	free(z_geom); 
        for(i=0;i<num_entries;++i)
	    free(felement[i]);
        free(simples);
	return; }

    void compute_B(void);
    void cart_to_internal(double*);
    void print_internals(void);
    void initial_H(void);
    void write_file30();
    void read_file11();
    void opt_step();   
};






