/*! \file stretch.h
    \ingroup (OPTKING)
    \brief Enter brief description of file here 
*/

class stretch_class {
    int id;
    int A;
    int B;
    double value; /* length of bond */
    double *s_A;  /* The s vector for atom A (Xa-Xb) */
    double *s_B;  /* The s vector for atom B (Xb-Xa) */
  public:
    stretch_class(){
      s_A = new double[3];
      s_B = new double[3];
    }
    ~stretch_class() {
      // fprintf(stdout,"destructing stretch class\n");
      delete [] s_A ;
      delete [] s_B ;
    }
    void print(FILE *fp_out, int print_flag) {
      if (print_flag == 0)
        fprintf(fp_out,"    (%d %d %d)\n", id, A+1, B+1);
      else
        fprintf(fp_out,"    (%d %d %d) (%.8lf)\n", id,A+1,B+1,value);
    }
    void set_id(int i) { id = i;}
    int  get_id(void)  { return id;}
    void set_A(int i)  { A = i;}
    int  get_A(void)   { return A;}
    void set_B(int i)  { B = i;}
    int  get_B(void)   { return B;}
    void    set_value(double length) { value = length;}
    double  get_value(void)          { return value;}
    void set_s_A(double s_A0, double s_A1, double s_A2) {
         s_A[0] = s_A0; s_A[1] = s_A1; s_A[2] = s_A2; }
    double get_s_A(int i) { return s_A[i]; }
    void set_s_B(double s_B0, double s_B1, double s_B2) {
         s_B[0] = s_B0; s_B[1] = s_B1; s_B[2] = s_B2; }
    double get_s_B(int i) { return s_B[i]; }
};




/*** STRETCH_SET class declaration ***/ 

class stretch_set {

  int num;
   stretch_class* stre_array;

  public:

   stretch_set(int size){
	if (0 <= size <10000){
	    stre_array = new stretch_class[size];
	  }
	else { fprintf(outfile,"\nWARNING: bad number of stretches\n"); }
      }

   stretch_set(void) { } /* don't allocate memory yet */
   void allocate(int size) {
	 if (0 <= size <10000)
	   stre_array = new stretch_class[size];
	 else 
       fprintf(outfile,"\nWARNING: bad number of stretches\n");
   }

   ~stretch_set() {
     // ~stretch_class is called automatically
     // fprintf(stdout,"destructing stretch_set\n");
     delete [] stre_array;
   }

   void print(FILE *fp_out, int print_flag) {
      int i;
      if (num > 0) {
        if (print_flag == 0) fprintf(fp_out,"  stre = (\n");
        else fprintf(fp_out,"Stretches\n");
        for (i=0; i < num; ++i)
           stre_array[i].print(fp_out, print_flag);
        if (print_flag == 0) fprintf(fp_out,"  )\n");
      }
      return;
    }
    void print_s() {
      int i;
      for (i=0;i<num;++i) {
        fprintf(outfile,"S vector for stretch %d %d: atom A\n",
         get_A(i),get_B(i) );
        fprintf(outfile,"(%16.10f,%16.10f,%16.10f)\n",
         get_s_A(i,0), get_s_A(i,1), get_s_A(i,2) );
        fprintf(outfile,"S vector for stretch %d %d: atom B\n",
         get_A(i),get_B(i) );
        fprintf(outfile,"(%16.10f,%16.10f,%16.10f)\n",
         get_s_B(i,0), get_s_B(i,1), get_s_B(i,2) );
      }
      return;
    }
    void compute(int natom, double *geom) {
      int i,A,B;
      double tmp;
      for (i=0; i< num ; ++i) {
        A = get_A(i);
        B = get_B(i);
        tmp = SQR( geom[3*A+0] - geom[3*B+0] ) +
              SQR( geom[3*A+1] - geom[3*B+1] ) +
              SQR( geom[3*A+2] - geom[3*B+2] );

        set_val(i, sqrt(tmp)*_bohr2angstroms);
      }
      return;
    }
    void compute_s(int natom, double *geom) {
      int i,j,A,B;
      double eBA[3], norm;
      for (i=0;i<num;++i) {
        A = get_A(i);
        B = get_B(i);

        for (j=0;j<3;++j)
          eBA[j] = geom[3*A+j] - geom[3*B+j];

        norm = sqrt( SQR(eBA[0]) + SQR(eBA[1]) + SQR(eBA[2]) );
        scalar_div(norm,eBA);

        set_s_A(i,eBA[0],eBA[1],eBA[2]);
        set_s_B(i, -1.0*eBA[0], -1.0*eBA[1], -1.0*eBA[2]);
      }
      return;
    }
    void set_num(int i) { num = i;}
    int  get_num(void) { return num;}
    void set_id(int index, int new_id) { stre_array[index].set_id(new_id);}
    int  get_id(int index) { return stre_array[index].get_id();}
    void set_A(int index, int new_A) { stre_array[index].set_A(new_A);}
    int  get_A(int index) {return stre_array[index].get_A();}
    void set_B(int index, int new_B) { stre_array[index].set_B(new_B);}
    int  get_B(int index) { return stre_array[index].get_B();}
    void set_val(int index, double new_val) { stre_array[index].set_value(new_val);}
    double  get_val(int index) { return stre_array[index].get_value();}

    void set_s_A(int index, double s_A0, double s_A1, double s_A2) {
                 stre_array[index].set_s_A(s_A0,s_A1,s_A2); }
    double get_s_A(int index, int i) { return stre_array[index].get_s_A(i); }

    void set_s_B(int index, double s_B0, double s_B1, double s_B2) {
                 stre_array[index].set_s_B(s_B0,s_B1,s_B2); }
    double get_s_B(int index, int i) { return stre_array[index].get_s_B(i); }
    int get_id_from_atoms(int a, int b) {
       int i;
  //     fprintf(outfile,"stre.get_id_from_atoms(%d,%d)\n",a,b);
       for (i=0;i<num;++i) {
         if ( (a == get_A(i)) && (b == get_B(i)) ) break;
       }
       if (i == num) {
         fprintf(outfile,"Could not find simple internal for atoms  \
             %d %d in list.\n", a+1, b+1);
         exit(2);
       }
   //    fprintf(outfile,"returning id: %d\n", get_id(i));
       return get_id(i);
    }
};

