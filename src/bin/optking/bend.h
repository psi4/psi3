
class bend_class {
    int id;
    int A;
    int B;
    int C;
    double value;
    double *s_A; // The s vector for atom A
    double *s_B;
    double *s_C;
  public:
    bend_class() {
      s_A = new double[3];
      s_B = new double[3];
      s_C = new double[3];
    }
    ~bend_class() {
    //  fprintf(stdout,"destructing bend class\n");
      delete [] s_A;
      delete [] s_B;
      delete [] s_C;
    }
    void print(FILE *fp_out, int print_flag) {
      if (print_flag == 0)
        fprintf(fp_out,"    (%d %d %d %d)\n", id,A+1,B+1,C+1);
      else 
        fprintf(fp_out,"    (%d %d %d %d) (%.8lf)\n", id,A+1,B+1,C+1,value);
    }
    void set_id(int i){ id = i;}
    int  get_id(void) { return id;}
    void set_A(int i) { A = i;}
    int  get_A(void)  { return A;}
    void set_B(int i) { B = i;}
    int  get_B(void)  { return B;}
    void set_C(int i) { C = i;}
    int  get_C(void)  { return C;}
    void set_value(double new_val) { value = new_val;}
    double get_value(void)  { return value;}
    void set_s_A(double s_A0, double s_A1, double s_A2) {
         s_A[0] = s_A0; s_A[1] = s_A1; s_A[2] = s_A2; }
    double get_s_A(int i) { return s_A[i]; }
    void set_s_B(double s_B0, double s_B1, double s_B2) {
         s_B[0] = s_B0; s_B[1] = s_B1; s_B[2] = s_B2; }
    double get_s_B(int i) { return s_B[i]; }
    void set_s_C(double s_C0, double s_C1, double s_C2) {
         s_C[0] = s_C0; s_C[1] = s_C1; s_C[2] = s_C2; }
    double get_s_C(int i) { return s_C[i]; }
};


class bend_set {

   int num;
   bend_class *bend_array;

  public:

   bend_set(int size) {
       if(0 <= size < 10000) {
	   bend_array = new bend_class[size];
	 }
       else { fprintf(outfile,"\nWARNING: bad number of bends\n");}
     }

   ~bend_set() {
     // fprintf(stdout,"destructing bend_set\n");
     delete [] bend_array;
   }

   void print(FILE *fp_out, int print_flag) {
      int i;
      if (num > 0) {
        if (print_flag == 0) fprintf(fp_out,"  bend = (\n");
        else fprintf(fp_out, "Bends\n");
        for (i=0; i < num; ++i)
           bend_array[i].print(fp_out, print_flag);
        if (print_flag == 0) fprintf(fp_out,"  )\n");
      }
      return;
    }

   void print_s() {
      int i;
      for (i=0;i<num;++i) {
        fprintf(outfile,"S vector for bend %d %d %d: atom A\n",get_A(i),get_B(i),get_C(i) );
        fprintf(outfile,"(%16.10f,%16.10f,%16.10f)\n", get_s_A(i,0), get_s_A(i,1), get_s_A(i,2) );
        fprintf(outfile,"S vector for bend %d %d %d: atom B\n",get_A(i),get_B(i),get_C(i) );
        fprintf(outfile,"(%16.10f,%16.10f,%16.10f)\n", get_s_B(i,0), get_s_B(i,1), get_s_B(i,2) );
        fprintf(outfile,"S vector for bend %d %d %d: atom C\n",get_A(i),get_B(i),get_C(i) );
        fprintf(outfile,"(%16.10f,%16.10f,%16.10f)\n", get_s_C(i,0), get_s_C(i,1), get_s_C(i,2) );
      }
      return;
    }

    void set_num(int i) { num = i;}
    int  get_num(void) { return num;}
    void set_id(int index, int new_id) { bend_array[index].set_id(new_id);}
    int  get_id(int index) { return bend_array[index].get_id();}
    void set_A(int index, int new_A) { bend_array[index].set_A(new_A);}
    int  get_A(int index) {return bend_array[index].get_A();}
    void set_B(int index, int new_B) { bend_array[index].set_B(new_B);}
    int  get_B(int index) { return bend_array[index].get_B();}
    void set_C(int index, int new_C) { bend_array[index].set_C(new_C);}
    int  get_C(int index) { return bend_array[index].get_C();}
    void set_val(int index, double new_val) { bend_array[index].set_value(new_val);}
    double  get_val(int index) { return bend_array[index].get_value();}

    void set_s_A(int index, double s_A0, double s_A1, double s_A2) {
                 bend_array[index].set_s_A(s_A0,s_A1,s_A2); }
    double get_s_A(int index, int i) { return bend_array[index].get_s_A(i); }

    void set_s_B(int index, double s_B0, double s_B1, double s_B2) {
                 bend_array[index].set_s_B(s_B0,s_B1,s_B2); }
    double get_s_B(int index, int i) { return bend_array[index].get_s_B(i); }

    void set_s_C(int index, double s_C0, double s_C1, double s_C2) {
                 bend_array[index].set_s_C(s_C0,s_C1,s_C2); }
    double get_s_C(int index, int i) { return bend_array[index].get_s_C(i); }
    void compute(int natom, double *geom) {
      int i,j,A,B,C;
      double rBA,rBC,eBA[3],eBC[3],tmp[3],dotprod,angle;
    
      for (i=0;i<num;++i) {
        A = get_A(i);
        B = get_B(i);
        C = get_C(i);
    
        for (j=0;j<3;++j) {
          eBA[j] = geom[3*A+j] - geom[3*B+j];
          eBC[j] = geom[3*C+j] - geom[3*B+j];
        }
    
        rBA = sqrt( SQR(eBA[0])+SQR(eBA[1])+SQR(eBA[2]) );
        rBC = sqrt( SQR(eBC[0])+SQR(eBC[1])+SQR(eBC[2]) );
    
        scalar_div(rBA,eBA);
        scalar_div(rBC,eBC);
    
        dot_arr(eBA,eBC,3,&dotprod);
    
        if (dotprod > 1.0) angle = 0.0;
        else if (dotprod < -1.0) angle = _pi;
        else angle = acos(dotprod)*180.0/_pi;
    
        set_val(i,angle);
      }
      return;
    }
    void compute_s(int natom, double *geom) {
      int i,j,A,B,C;
      double val,rBA,rBC;
      double eBA[3], eBC[3], tmp[3];
      double *geom_ang;
    
      geom_ang  = new double[3*natom];
      for (i=0;i<natom*3;++i)
        geom_ang[i] = geom[i] * _bohr2angstroms;
    
      for (i=0;i<num;++i) {
        A = get_A(i);
        B = get_B(i);
        C = get_C(i);
        val = get_val(i)*_pi/180.0;
    
        for (j=0;j<3;++j) {
          eBA[j] = geom_ang[3*A+j] - geom_ang[3*B+j];
          eBC[j] = geom_ang[3*C+j] - geom_ang[3*B+j];
        }
    
        rBA = sqrt( SQR(eBA[0]) + SQR(eBA[1]) + SQR(eBA[2]) );
        rBC = sqrt( SQR(eBC[0]) + SQR(eBC[1]) + SQR(eBC[2]) );

        for (j=0;j<3;++j) {
          eBA[j] = eBA[j] / rBA;
          eBC[j] = eBC[j] / rBC;
        }

        for (j=0;j<3;++j) {
          tmp[j] = (eBA[j]*cos(val) - eBC[j]) / (rBA*sin(val));
        }
        set_s_A(i,tmp[0],tmp[1],tmp[2]);
    
        for (j=0;j<3;++j) {
          tmp[j] = ((rBA - rBC*cos(val))*eBA[j] + (rBC-rBA*cos(val))*eBC[j])
                    / (rBA * rBC * sin(val));
        }
        set_s_B(i,tmp[0],tmp[1],tmp[2]);
    
        for (j=0;j<3;++j) {
          tmp[j] = (eBC[j]*cos(val) - eBA[j]) / (rBC*sin(val));
        }
        set_s_C(i,tmp[0],tmp[1],tmp[2]);
      }
      delete [] geom_ang;
      return;
    }
    int get_id_from_atoms(int a, int b, int c) {
       int i;
//       fprintf(outfile,"bend.get_id_from_atoms(%d,%d,%d)\n",a,b,c);
       for (i=0;i<num;++i) {
         if ( (a == get_A(i)) && (b == get_B(i)) && (c == get_C(i)) ) break;
       }
 //      fprintf(outfile,"Returning id: %d\n",get_id(i));
       return get_id(i);
    }
};

