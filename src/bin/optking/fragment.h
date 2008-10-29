/*! \file
    \ingroup OPTKING
    \brief Six coordinates for non-bonded fragments
*/

#ifndef _psi3_bin_optking_fragment_h_
#define _psi3_bin_optking_fragment_h_

namespace psi { namespace optking {

class fragment_class {
    int id;
    int J;         /* which of six fragment types */
    int A_natom;   /* number of atoms in fragment A */
    int B_natom;   /* number of atoms in fragment B */
    int *A_atom;   /* list of atoms in fragment A */
    int *B_atom;   /* list of atoms in fragment B */
    double value;  /* D_J = value of coordinate */
    double *A_weight; /* weights can be used to fix 3 points using a fixed linear */
    double *B_weight; /* combination of atom locations, dimension is: 3 x A_natom */
    double *A_s;  /* S vectors for J in fragment A */
    double *B_s;  /* S vectors for J in fragment B */
//    int P;        /* the number of reference points for the fragment - 
//                     P=1 (for an atom); P=2 (for a linear fragment); P=3 (for non-linear fragments) */
//    int dim;      /* the number of interfragment coordinates -
//                     dim=1 (for an atom); 

  public:

    fragment_class() { }

    void allocate() {
      A_atom   = new int[A_natom];
      B_atom   = new int[B_natom];
      // A_weight and B_weight are allocated in internals
      A_s = init_array(9); /* 3 reference points with xyz */
      B_s = init_array(9);
    }

    ~fragment_class() {
      delete [] A_atom;
      delete [] B_atom;
      delete [] A_weight;
      delete [] B_weight;
      free(A_s);
      free(B_s);
    }

    /* functions in fragment.cc */
    void print(FILE *fp_out, int print_flag, int print_weights);
    void compute(double *geom);
    void compute_s(int natom, double *geom);
    double get_val_A_or_rad(void);
    void print_s(void);

    /* simple functions */
    void set_id(int i){ id = i;}
    int  get_id(void) { return id;}
    void set_J(int j){ J = j;}
    int  get_J(void) { return J;}
    void set_A_natom(int i) { A_natom = i;}
    int  get_A_natom(void)  { return A_natom;}
    void set_B_natom(int i) { B_natom = i;}
    int  get_B_natom(void)  { return B_natom;}
    void   set_value(double new_value) { value = new_value;}
    double get_value(void)  { return value;}
    void  set_A_atom(int i, int j) { A_atom[i] = j;}
    int   get_A_atom(int i) { return A_atom[i];}
    void  set_B_atom(int i, int j) { B_atom[i] = j;}
    int   get_B_atom(int i) { return B_atom[i];}

    void   set_A_weight(double *new_A_weight) { A_weight = new_A_weight; }
    double get_A_weight(int ref_atom, int frag_index) { return A_weight[A_natom*ref_atom+frag_index]; }
    void   set_B_weight(double *new_B_weight) { B_weight = new_B_weight; }
    double get_B_weight(int ref_atom, int frag_index) { return B_weight[B_natom*ref_atom+frag_index]; }

    double get_A_s(int ref_atom_xyz) { return A_s[ref_atom_xyz]; }
    double get_B_s(int ref_atom_xyz) { return B_s[ref_atom_xyz]; }
};

class fragment_set {

   int num;
   fragment_class *fragment_array;

  public:

   fragment_set(int size) {
       if(0 <= size < 1000) {
	   fragment_array = new fragment_class[size];
	 }
       else { fprintf(outfile,"\nWARNING: bad number of fragments\n");}
     }

   fragment_set(void) { } /* don't allocate memory yet */

   void allocate(int size) {
     int i;
     if (0 <= size <10000)
       fragment_array = new fragment_class[size];
     else
       fprintf(outfile,"\nWARNING: bad number of fragments\n");
   }

   ~fragment_set() {
     delete [] fragment_array;
   }

   void allocate_one(int index) { fragment_array[index].allocate(); }

   void print(FILE *fp_out, int print_flag, int print_weights) {
      int i;
      if (num > 0) {
        if (print_flag == 0) fprintf(fp_out,"  fragment = (\n");
        else fprintf(fp_out, "Fragments\n");
        for (i=0; i < num; ++i)
           fragment_array[i].print(fp_out, print_flag, print_weights);
        if (print_flag == 0) fprintf(fp_out,"  )\n");
      }
      return;
    }

    void compute(double *geom) {
      int i;
      for (i=0;i<num;++i)
        fragment_array[i].compute(geom);
      return;
    }

    void compute_s(int natom, double *geom) {
      int i;
      for (i=0;i<num;++i)
        fragment_array[i].compute_s(natom, geom);
      return;
    }
    void print_s() {
      int i;
      for (i=0;i<num;++i)
        fragment_array[i].print_s();
    }

    void set_num(int i) { num = i;}
    int  get_num(void) { return num; } // may want a more elegant solution later
    void set_id(int index, int new_id) { fragment_array[index].set_id(new_id);}
    int  get_id(int index) { return fragment_array[index].get_id();}
    void set_J(int index, int new_J) { fragment_array[index].set_J(new_J);}
    int  get_J(int index) { return fragment_array[index].get_J();}
    void set_A_natom(int index, int new_A_natom) { fragment_array[index].set_A_natom(new_A_natom);}
    int  get_A_natom(int index) {return fragment_array[index].get_A_natom();}
    void set_B_natom(int index, int new_B_natom) { fragment_array[index].set_B_natom(new_B_natom);}
    int  get_B_natom(int index) {return fragment_array[index].get_B_natom();}
    void   set_val(int index, double new_val) { fragment_array[index].set_value(new_val);}
    double get_val(int index) { return fragment_array[index].get_value();}
    double get_val_A_or_rad(int index) { return fragment_array[index].get_val_A_or_rad();}

    void set_A_atom(int index, int frag_index, int atom)
      { fragment_array[index].set_A_atom(frag_index, atom);}
    int  get_A_atom(int index, int frag_index)
      { return fragment_array[index].get_A_atom(frag_index);}
    void set_B_atom(int index, int frag_index, int atom)
      { fragment_array[index].set_B_atom(frag_index, atom);}
    int  get_B_atom(int index, int frag_index)
      { return fragment_array[index].get_B_atom(frag_index);}

    void   set_A_weight(int index, double *new_A_weight)
      { fragment_array[index].set_A_weight(new_A_weight); }
    double get_A_weight(int index, int ref_atom, int frag_index)
      { return fragment_array[index].get_A_weight(ref_atom, frag_index); }
    void   set_B_weight(int index, double *new_B_weight)
      { fragment_array[index].set_B_weight(new_B_weight); }
    double get_B_weight(int index, int ref_atom, int frag_index)
      { return fragment_array[index].get_B_weight(ref_atom, frag_index); }

    double get_A_s(int index, int atom_xyz)
      { return fragment_array[index].get_A_s(atom_xyz); }
    double get_B_s(int index, int atom_xyz)
      { return fragment_array[index].get_B_s(atom_xyz); }

    int get_id_from_atoms(int a_natom, int b_natom, int *a_atom, int *b_atom, int j) {
      int i, a, b, match=0;
      while ((match==0) && (i<num)) {
         match = 1;
         if ( (a_natom == get_A_natom(i)) && (b_natom == get_B_natom(i)) ) {
           for (a=0; a<a_natom; ++a) {
             if (a_atom[a] != get_A_atom(i,a))
               match=0;
           }
           for (b=0; b<b_natom; ++b) {
             if (b_atom[b] != get_B_atom(i,b))
               match=0;
           }
           if (get_J(i) != j) match = 0;
         }
         else match = 0;
         ++i;
       }
       
       if ((i == num) || (match == 0)) {
         fprintf(outfile,"Could not find simple fragment with natoms %d %d and J=%d.\n",
          a_natom, b_natom, j);
         exit(2);
       }
       return get_id(i-1);
    }
};

}} /* namespace psi::optking */

#endif
