/*##################################################################################
  simple_internal.h

  simple coordinate concrete class declaration and definition
  a type which coord_base can hold an array of

  we need a consistent definition for the internal type:
     bond length = 0
     dihedral angle = 1;
     torsional angle = 2;
  ################################################################################*/

class simple_internal {
    int type;
    double val;
  public:
    void set(int i) { type = i; return;}
    void set(double a) { val = a; return;}
    void set(int i, double a) { type = i; val = a; return;}
    void set(double a, int i) { type = i; val = a; return;}
    int get_type() { return type;}
    double get_val() { return val;}
};
    
    
