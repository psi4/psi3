/* class for pairs of orthogonal linear bend coordinates
 * lin with linval1 is toward an axis perpendicular with A-B-C with maximum
 * overlap with the y-axis
 * lin with linval2 is orthogonal to this axis
 */

extern "C" {
  #include <libchkpt/chkpt.h>
}

class lin_bend_class {
    int id;
    int A;
    int B;
    int C;
    int linval; // 1 or 2 for lin1 or lin2
    double value;
    double *s_A; /* The s vector for atom A */
    double *s_B;
    double *s_C;
    double *dummy; // x,y,z of implicit dummy atom to orient lin_bend
  public:
    lin_bend_class() {
      s_A = new double[3];
      s_B = new double[3];
      s_C = new double[3];
      dummy = new double[3];
    }
    ~lin_bend_class() {
    // fprintf(stdout,"destructing lin_bend class\n"); fflush(outfile);
      delete [] s_A;
      delete [] s_B;
      delete [] s_C;
      delete [] dummy;
    }
    void print(FILE *fp_out, int print_flag) {
      if (print_flag == 0) // to intco.dat - don't give linval
        fprintf(fp_out,"    (%d %d %d %d)\n", id,A+1,B+1,C+1);
      else 
        fprintf(fp_out,"    (%d %d %d %d %d) (%.8lf)\n", id,A+1,B+1,C+1,linval,value);
    }
    void set_id(int i){ id = i;}
    int  get_id(void) { return id;}
    void set_A(int i) { A = i;}
    int  get_A(void)  { return A;}
    void set_B(int i) { B = i;}
    int  get_B(void)  { return B;}
    void set_C(int i) { C = i;}
    int  get_C(void)  { return C;}
    void set_linval(int i){ linval = i;}
    int  get_linval(void) { return linval;}
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
    double get_dummy(int i) { return dummy[i]; }
    void set_dummy(double x, double y, double z) {
         dummy[0] = x;  dummy[1] = y; dummy[2] = z; }
};


class lin_bend_set {

  int num; // lin1,lin2 pair count as 2
  lin_bend_class *lin_bend_array;

  public:

  lin_bend_set(int size) {
    if(0 <= size < 100)
      lin_bend_array = new lin_bend_class[size];
    else
      fprintf(outfile,"\nWARNING: bad number of linear bends.\n");
  }

  lin_bend_set(void) { } /* don't allocate memory yet */
  void allocate(int size) {
    if (0 <= size <10000)
      lin_bend_array = new lin_bend_class[size];
    else
      fprintf(outfile,"\nWARNING: bad number of lin_bends\n");
  }

  ~lin_bend_set() {
    /* fprintf(stdout,"destructing lin_bend_set\n"); */
    delete [] lin_bend_array;
  }

  // print_flag = 0 goes to intco.dat
  // print_flag = 1 goes to output.dat
  void print(FILE *fp_out, int print_flag) {
    int i;
    if (num > 0) {
      if (print_flag == 0) {
        fprintf(fp_out,"  lin1 = (\n");
        for (i=0; i < num; ++i) {
          if (get_linval(i) == 1)
            lin_bend_array[i].print(fp_out, print_flag);
        }
        fprintf(fp_out,"  )\n");
        fprintf(fp_out,"  lin2 = (\n");
        for (i=0; i < num; ++i) {
          if (get_linval(i) == 2)
            lin_bend_array[i].print(fp_out, print_flag);
        }
        fprintf(fp_out,"  )\n");
      }
      else {
        fprintf(fp_out, "Linear Bends\n");
        for (i=0; i < num; ++i)
          lin_bend_array[i].print(fp_out, print_flag);
      }
    }
    return;
  }

  void print_s() {
    int i;
    for (i=0;i<num;++i) {
      fprintf(outfile,"S vector for lin_bend %d %d %d %d: atom A\n",
          get_A(i),get_B(i),get_C(i),get_linval(i) );
      fprintf(outfile,"(%16.10f,%16.10f,%16.10f)\n",
          get_s_A(i,0), get_s_A(i,1), get_s_A(i,2) );
      fprintf(outfile,"S vector for lin_bend %d %d %d %d: atom B\n",
          get_A(i),get_B(i),get_C(i),get_linval(i) );
      fprintf(outfile,"(%16.10f,%16.10f,%16.10f)\n",
          get_s_B(i,0), get_s_B(i,1), get_s_B(i,2) );
      fprintf(outfile,"S vector for lin_bend %d %d %d %d: atom C\n",
          get_A(i),get_B(i),get_C(i),get_linval(i) );
      fprintf(outfile,"(%16.10f,%16.10f,%16.10f)\n",
          get_s_C(i,0), get_s_C(i,1), get_s_C(i,2) );
    }
    return;
  }

  void set_num(int i) { num = i;}
  int  get_num(void) { return num;}
  void set_id(int index, int new_id) { lin_bend_array[index].set_id(new_id);}
  int  get_id(int index) { return lin_bend_array[index].get_id();}
  void set_A(int index, int new_A) { lin_bend_array[index].set_A(new_A);}
  int  get_A(int index) {return lin_bend_array[index].get_A();}
  void set_B(int index, int new_B) { lin_bend_array[index].set_B(new_B);}
  int  get_B(int index) { return lin_bend_array[index].get_B();}
  void set_C(int index, int new_C) { lin_bend_array[index].set_C(new_C);}
  int  get_C(int index) { return lin_bend_array[index].get_C();}
  void set_linval(int index, int new_linval) { lin_bend_array[index].set_linval(new_linval);}
  int  get_linval(int index) { return lin_bend_array[index].get_linval();}
  void set_val(int index, double new_val) { lin_bend_array[index].set_value(new_val);}
  double  get_val(int index) { return lin_bend_array[index].get_value();}
  void set_s_A(int index, double s_A0, double s_A1, double s_A2) {
               lin_bend_array[index].set_s_A(s_A0,s_A1,s_A2); }
  double get_s_A(int index, int i) { return lin_bend_array[index].get_s_A(i); }

  void set_s_B(int index, double s_B0, double s_B1, double s_B2) {
               lin_bend_array[index].set_s_B(s_B0,s_B1,s_B2); }
  double get_s_B(int index, int i) { return lin_bend_array[index].get_s_B(i); }

  void set_s_C(int index, double s_C0, double s_C1, double s_C2) {
               lin_bend_array[index].set_s_C(s_C0,s_C1,s_C2); }
  double get_s_C(int index, int i) { return lin_bend_array[index].get_s_C(i); }

  void set_dummy(int index, double d_x, double d_y, double d_z) {
               lin_bend_array[index].set_dummy(d_x,d_y,d_z); }
  double get_dummy(int index, int i) { return lin_bend_array[index].get_dummy(i); }


  void compute(int natom, double *geom) {
    int i,j,A,B,C,linval,rottype;
    double rBA,rBC,rBD,eBA[3],eBC[3],eBD[3],tmp[3],dotprod;
    double dummy[3], angle_ABD, angle_CBD;
    
    for (i=0;i<num;++i) {
      A = get_A(i);
      B = get_B(i);
      C = get_C(i);
      linval = get_linval(i);

      chkpt_init(PSIO_OPEN_OLD);
      rottype = chkpt_rd_rottype();
      chkpt_close();

      // if we zoom dummy atom way out, then our s vectors can just
      // be steps toward or away from the dummy atom
      // the following implicit dummy atom specification is
      // certainly not a general solution
     
      if (rottype == 0) { // assymmetric top - use z-axis
        set_dummy(i, 0.0, 0.0, 1.0E6 );
      }
      else { // otherwise use x-axis; not a general solution
        set_dummy(i, 1.0E6, 0.0, 0.0 );
      }

      dummy[0] = get_dummy(i,0);
      dummy[1] = get_dummy(i,1);
      dummy[2] = get_dummy(i,2);

      // the following will compute the placement of the other
      // implicit dummy atom in general - may use later
      if (linval == 2) {
        //positive displacement is toward eBD X eBA
        for (j=0;j<3;++j) {
          eBD[j] = dummy[j] - geom[3*B+j];
          eBA[j] = geom[3*A+j] - geom[3*B+j];
        }
        rBD = sqrt( SQR(eBD[0])+SQR(eBD[1])+SQR(eBD[2]) );
        rBA = sqrt( SQR(eBA[0])+SQR(eBA[1])+SQR(eBA[2]) );
        scalar_div(rBD,eBD);
        scalar_div(rBA,eBA);
        cross_product(eBD,eBA,tmp);

        set_dummy(i,1.0E6*tmp[0],1.0E6*tmp[1],1.0E6*tmp[2]);
        dummy[0] = get_dummy(i,0);
        dummy[1] = get_dummy(i,1);
        dummy[2] = get_dummy(i,2);
      }

      //fprintf(outfile,"dummy atom at %10.5lf %10.5lf %10.5lf\n",
          //dummy[0],dummy[1],dummy[2]);

      // angle = <ABD + <CBD
      // compute value of A-B-D
      for (j=0;j<3;++j) {
        eBA[j] = geom[3*A+j] - geom[3*B+j];
        eBC[j] = dummy[j] - geom[3*B+j];
      }
      rBA = sqrt( SQR(eBA[0])+SQR(eBA[1])+SQR(eBA[2]) );
      rBC = sqrt( SQR(eBC[0])+SQR(eBC[1])+SQR(eBC[2]) );
      scalar_div(rBA,eBA);
      scalar_div(rBC,eBC);
      dot_arr(eBA,eBC,3,&dotprod);
      if (dotprod > 1.0) angle_ABD = 0.0;
      else if (dotprod < -1.0) angle_ABD = _pi;
      else angle_ABD = acos(dotprod)*180.0/_pi;

      // compute value of CBD
      for (j=0;j<3;++j) {
        eBA[j] = geom[3*C+j] - geom[3*B+j];
        eBC[j] = dummy[j] - geom[3*B+j];
      }
      rBA = sqrt( SQR(eBA[0])+SQR(eBA[1])+SQR(eBA[2]) );
      rBC = sqrt( SQR(eBC[0])+SQR(eBC[1])+SQR(eBC[2]) );
      scalar_div(rBA,eBA);
      scalar_div(rBC,eBC);
      dot_arr(eBA,eBC,3,&dotprod);
      if (dotprod > 1.0) angle_CBD = 0.0;
      else if (dotprod < -1.0) angle_CBD = _pi;
      else angle_CBD = acos(dotprod)*180.0/_pi;

      set_val(i,angle_ABD+angle_CBD);
    } // end loop over linear angles
    return;
  }


  // s vectors point in direction of increasing internal coordinate value
  // so A and C retreat from D and B advances
  void compute_s(int natom, double *geom) {
    int i,j,A,B,C;
    double rAD, rBD, rCD, eAD[3], eBD[3], eCD[3], dummy[3];

    for (i=0;i<num;++i) {
      dummy[0] = get_dummy(i,0);
      dummy[1] = get_dummy(i,1);
      dummy[2] = get_dummy(i,2);

      A = get_A(i);
      B = get_B(i);
      C = get_C(i);

      // zoom D way out along BD 
      for (j=0;j<3;++j)
        eBD[j] = dummy[j] - geom[3*B+j];

      rBD = sqrt( SQR(eBD[0])+SQR(eBD[1])+SQR(eBD[2]) );
      scalar_div(rBD,eBD);
      for (j=0;j<3;++j)
        dummy[j] += eBD[j] * 1.0E9; 

      for (j=0;j<3;++j) {
        eAD[j] = dummy[j] - geom[3*A+j];
        eBD[j] = dummy[j] - geom[3*B+j];
        eCD[j] = dummy[j] - geom[3*C+j];
      }
      rAD = sqrt( SQR(eAD[0])+SQR(eAD[1])+SQR(eAD[2]) );
      rBD = sqrt( SQR(eBD[0])+SQR(eBD[1])+SQR(eBD[2]) );
      rCD = sqrt( SQR(eCD[0])+SQR(eCD[1])+SQR(eCD[2]) );
      scalar_div(rAD,eAD);
      scalar_div(rBD,eBD);
      scalar_div(rCD,eCD);
  
      // A and C go away from D
      set_s_A(i, -eAD[0], -eAD[1], -eAD[2]);
      set_s_C(i, -eCD[0], -eCD[1], -eCD[2]);
  
      // B goes toward D
      set_s_B(i, eBD[0], eBD[1], eBD[2]);
    }
    return;
  }

  int get_id_from_atoms(int a, int b, int c, int linval) {
    int i;
      /* fprintf(outfile,"lin_bend.get_id_from_atoms(%d,%d,%d)\n",a,b,c); */
     for (i=0;i<num;++i) {
       if ( (a == get_A(i)) && (b == get_B(i)) && (c == get_C(i))
         && (linval == get_linval(i)) ) break;
     }
     /* fprintf(outfile,"Returning id: %d\n",get_id(i)); */
     return get_id(i);
  }
};

