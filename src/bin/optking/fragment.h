/*! \file
    \ingroup OPTKING
    \brief fragment_class describes a set of interfragment coordinates connecting two fragments
*/

#ifndef _psi3_bin_optking_fragment_h_
#define _psi3_bin_optking_fragment_h_

namespace psi { namespace optking {

/*
Fragment_class describes a SET of interfragment coordinates connecting two fragments
(A and B).  Each member of the set shares the definition of the fragments:
(A_natom, *A_atom, *A_weight, B_natom, *B_atom, *B_weight) and id number.
Each member gets its own s-vectors.
The following index order of values and s-vectors (A_s and B_s) is assumed:
RAB, theta_A, theta_B, tau, phi_A, phi_B.
*/

class fragment_class {
    bool *coord_on; /* indicates on/off for which coordinates [1-6] are present */
    int id;        /* unique id number for the set */
    int A_natom;   /* number of atoms in fragment A */
    int B_natom;   /* number of atoms in fragment B */
    int *A_atom;   /* list of atoms in fragment A */
    int *B_atom;   /* list of atoms in fragment B */
    /* P=1 (for an atom); P=2 (for a linear fragment); P=3 (for non-linear fragments) */
    int A_P;       /* the # of reference points for fragment A to worry about */
    int B_P;       /* the # of reference points for fragment B to worry about; A_P >= B_P */
    double **A_weight; /* weights fix reference points via a linear combination [3][A_natom] */
    double **B_weight; /* weights fix reference points via a linear combination [3][B_natom] */
    double *value;  /* D = value of coordinate [6] */
    double **A_s; /* S vectors fragment A [6][A_P*3] */
    double **B_s; /* S vectors fragment B [6][B_P*3] */
    int *near_180; /* [6] ; +1 => val>160, -1 => value<160, 0 otherwise */

  public:

    fragment_class() { }

    void allocate() {
      A_atom   = new int[A_natom];
      B_atom   = new int[B_natom];
      coord_on = new bool[6];
      value    = new double[6];
      A_weight = block_matrix(3,A_natom);
      B_weight = block_matrix(3,B_natom);
      A_s = block_matrix(6,A_P*3);
      B_s = block_matrix(6,B_P*3);
      near_180 = new int[6];
    }

    ~fragment_class() {
      delete [] A_atom;
      delete [] B_atom;
      delete [] coord_on;
      delete [] value;
      free_block(A_weight);
      free_block(B_weight);
      free_block(A_s);
      free_block(B_s);
      delete [] near_180;
    }

    /* functions in fragment.cc */
    void print(FILE *fp_out, int print_flag, int print_weights);
    void compute(double *geom);
    void compute_s(int natom, double *geom);
    double get_val_A_or_rad(int I);
    void print_s(void);
    void fix_near_180(void); // in fragment.cc

    /*** set/get shared variables ***/
    void set_id(int i){ id = i;}
    int  get_id(void) { return id;}

    int  get_dim(void) {
      int I, dim=0;
      for (I=0; I<6; ++I) 
        if (coord_on[I]) ++dim;
      return dim;
    }

    void set_A_natom(int i) { A_natom = i;}
    int  get_A_natom(void)  { return A_natom;}
    void set_B_natom(int i) { B_natom = i;}
    int  get_B_natom(void)  { return B_natom;}

    void set_A_P(int i) { A_P = i;}
    int  get_A_P(void)  { return A_P;}
    void set_B_P(int i) { B_P = i;}
    int  get_B_P(void)  { return B_P;}

    void set_A_atom(int i, int j) { A_atom[i] = j;}
    int  get_A_atom(int i) {
      if (i >= A_natom) throw("fragment.get_A_atom() - index is too large");
      return A_atom[i];
    }
    void set_B_atom(int i, int j) { B_atom[i] = j;}
    int  get_B_atom(int i) {
      if (i >= B_natom) throw("fragment.get_B_atom() - index is too large");
      return B_atom[i];
    }

    void set_A_weight(int ref_atom, int frag_index, double new_weight) {
      if (ref_atom >= A_P)
        throw("fragment.get_A_weight() ref_atom is greater than A_P\n");
      else if (frag_index >= A_natom)
        throw("fragment.get_A_weight() frag_index is greater than A_natom\n");
      A_weight[ref_atom][frag_index] = new_weight;
    }
    void set_B_weight(int ref_atom, int frag_index, double new_weight) {
      if (ref_atom >= B_P)
        throw("fragment.get_B_weight() ref_atom is greater than B_P\n");
      else if (frag_index >= B_natom)
        throw("fragment.get_B_weight() frag_index is greater than B_natom\n");
      B_weight[ref_atom][frag_index] = new_weight;
    }

    /* ref_atom = 0-2, frag_index 0-3*A_natom */
    double get_A_weight(int ref_atom, int frag_index) {
      if (ref_atom >= A_P)
        throw("fragment.get_A_weight() ref_atom is greater than A_P\n");
      else if (frag_index >= A_natom)
        throw("fragment.get_A_weight() frag_index is greater than A_natom\n");
      return A_weight[ref_atom][frag_index];
    }
    double get_B_weight(int ref_atom, int frag_index) {
      if (ref_atom >= B_P)
        throw("fragment.get_B_weight() ref_atom is greater than B_P\n");
      else if (frag_index >= B_natom)
        throw("fragment.get_B_weight() frag_index is greater than B_natom\n");
      return B_weight[ref_atom][frag_index];
    }

    void set_near_180(int I, int new_val) {
      if (I<0 || I>5) throw("fragment.set_near_180() expects an id between 0 and 5");
      if ((new_val != -1) && (new_val != 1) && (new_val != 0))
        throw("fragment.set_near_180() does not understand new value");
      near_180[I] = new_val;
   }
   int get_near_180(int I) {
      if (I<0 || I>5) throw("fragment.get_near_180() expects an id between 0 and 5");
      return near_180[I];
   }

    /* set/get coordinate specific variables */
    void set_coord_on(int I, bool on_or_off) {
      if (I<0 || I>5) throw("fragment.set_coord_on() expects an id between 0 and 5");
      coord_on[I] = on_or_off;
    }
    bool get_coord_on(int I) {
      if (I<0 || I>5) throw("fragment.get_coord_on() expects an id between 0 and 5");
      return coord_on[I];
    }
    void set_value(int I, double new_value) {
      if (I<0 || I>5) throw("fragment.set_value() expects an id between 0 and 5");
      value[I] = new_value;
    }
    double get_value(int I)  {
      if (I<0 || I>5) throw("fragment.get_value() expects an id between 0 and 5");
      else if (!coord_on[I]) throw("fragment.get_value() - coordinate is not active");
      return value[I];
    }
    double get_A_s(int I, int ref_atom_xyz) {
      if (I<0 || I>5)
        throw("fragment.get_A_s() expects an id between 0 and 5");
      else if (ref_atom_xyz >= 3*A_P)
        throw("fragment.get_A_s() ref_atom_xyz >= 3*A_natom\n");
      return A_s[I][ref_atom_xyz];
    }
    double get_B_s(int I, int ref_atom_xyz) {
      if (I<0 || I>5)
        throw("fragment.get_B_s() expects an id between 0 and 5");
      else if (ref_atom_xyz >= 3*B_P)
        throw("fragment.get_B_s() ref_atom_xyz >= 3*B_natom\n");
      return B_s[I][ref_atom_xyz];
    }
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

    void fix_near_180(void) {
      int i;
      for (i=0; i<num; ++i)
        fragment_array[i].fix_near_180();
      return;
    }

    void print_s() {
      int i;
      for (i=0;i<num;++i)
        fragment_array[i].print_s();
    }

    void set_num(int i) {num = i;}
    int get_num(void) { return num;}

    int get_dim(void) {
      int i, dim=0;
      for (i=0; i<num; ++i)
        dim += fragment_array[i].get_dim();
      return dim;
    }

    int  get_dim(int index) { return fragment_array[index].get_dim(); }

    void set_id(int index, int new_id) { fragment_array[index].set_id(new_id);}
    int  get_id(int index) { return fragment_array[index].get_id();}

    void set_A_natom(int index, int new_A_natom) { fragment_array[index].set_A_natom(new_A_natom);}
    int  get_A_natom(int index) {return fragment_array[index].get_A_natom();}
    void set_B_natom(int index, int new_B_natom) { fragment_array[index].set_B_natom(new_B_natom);}
    int  get_B_natom(int index) {return fragment_array[index].get_B_natom();}

    void set_A_atom(int index, int frag_index, int atom)
      { fragment_array[index].set_A_atom(frag_index, atom);}
    int  get_A_atom(int index, int frag_index)
      { return fragment_array[index].get_A_atom(frag_index);}
    void set_B_atom(int index, int frag_index, int atom)
      { fragment_array[index].set_B_atom(frag_index, atom);}
    int  get_B_atom(int index, int frag_index)
      { return fragment_array[index].get_B_atom(frag_index);}

    void   set_A_weight(int index, int ref_atom, int frag_index, double new_weight)
      { fragment_array[index].set_A_weight(ref_atom, frag_index, new_weight); }
    void   set_B_weight(int index, int ref_atom, int frag_index, double new_weight)
      { fragment_array[index].set_B_weight(ref_atom, frag_index, new_weight); }
    double get_A_weight(int index, int ref_atom, int frag_index)
      { return fragment_array[index].get_A_weight(ref_atom, frag_index); }
    double get_B_weight(int index, int ref_atom, int frag_index)
      { return fragment_array[index].get_B_weight(ref_atom, frag_index); }

    int get_A_P(int index) { return fragment_array[index].get_A_P();}
    void set_A_P(int index, int i) { fragment_array[index].set_A_P(i);}
    int get_B_P(int index) { return fragment_array[index].get_B_P();}
    void set_B_P(int index, int i) { fragment_array[index].set_B_P(i);}

    /*** get/set I specific variables ***/
    void set_coord_on(int index, int I, bool on_or_off) {
      fragment_array[index].set_coord_on(I, on_or_off);
    }
    bool get_coord_on(int index, int I) {
      return fragment_array[index].get_coord_on(I);
    }

    void set_near_180(int index, int I, int new_val) {
      fragment_array[index].set_near_180(I, new_val);
    }
    int get_near_180(int index, int I) {
      return fragment_array[index].get_near_180(I);
    }

    void   set_val(int index, int I, double new_val) {
      fragment_array[index].set_value(I, new_val);
    }
    double get_val(int index, int I) {
      return fragment_array[index].get_value(I);
    }
    double get_val_A_or_rad(int index, int I) {
      return fragment_array[index].get_val_A_or_rad(I);
    }

    double get_A_s(int index, int I, int atom_xyz)
      { return fragment_array[index].get_A_s(I, atom_xyz); }
    double get_B_s(int index, int I, int atom_xyz)
      { return fragment_array[index].get_B_s(I, atom_xyz); }

    int get_id_from_atoms(int a_natom, int b_natom, int *a_atom, int *b_atom, int I) {
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
           if (get_coord_on(i,I) != I) match = 0;
         }
         else match = 0;
         ++i;
       }
       
       if ((i == num) || (match == 0)) {
         fprintf(outfile,"Could not find simple fragment with natoms %d %d and I=%d.\n",
          a_natom, b_natom, I);
         exit(2);
       }
       return get_id(i-1);
    }
    int *atom2fragment(int natom);

};

}} /* namespace psi::optking */

#endif
