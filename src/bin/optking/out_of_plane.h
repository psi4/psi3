class out_class {
    int id;
    int A;
    int B;
    int C;
    int D;
    double value;
    double s_A[3]; // The s vector for atom A
    double s_B[3];
    double s_C[3];
    double s_D[3];
  public:
    out_class(){}
    ~out_class(){}
    void print(FILE *fp_out, int print_flag) {
      if (print_flag == 0) 
        fprintf(fp_out,"    (%d %d %d %d %d)\n", id,A+1,B+1,C+1,D+1);
      else
        fprintf(fp_out,"    (%d %d %d %d %d) (%.8lf)\n", id,A+1,B+1,C+1,D+1,value);
    }
    void set_id(int i){ id = i;}
    int  get_id(void) { return id;}
    void set_A(int i) { A = i;}
    int  get_A(void)  { return A;}
    void set_B(int i) { B = i;}
    int  get_B(void)  { return B;}
    void set_C(int i) { C = i;}
    int  get_C(void)  { return C;}
    void set_D(int i) { D = i;}
    int  get_D(void)  { return D;}
    void set_value(double new_val) { value = new_val;}
    double get_value(void)  { return value;}
    void set_s_A(double s_A0, double s_A1, double s_A2) {
         s_A[0] = s_A0;
         s_A[1] = s_A1;
         s_A[2] = s_A2; }
    double get_s_A(int i) { return s_A[i]; }
    void set_s_B(double s_B0, double s_B1, double s_B2) {
         s_B[0] = s_B0;
         s_B[1] = s_B1;
         s_B[2] = s_B2; }
    double get_s_B(int i) { return s_B[i]; }
    void set_s_C(double s_C0, double s_C1, double s_C2) {
         s_C[0] = s_C0;
         s_C[1] = s_C1;
         s_C[2] = s_C2; }
    double get_s_C(int i) { return s_C[i]; }
    void set_s_D(double s_D0, double s_D1, double s_D2) {
         s_D[0] = s_D0;
         s_D[1] = s_D1;
         s_D[2] = s_D2; }
    double get_s_D(int i) { return s_D[i]; }
};

class out_set {

   int num;
   out_class *out_array;

  public:
 
   out_set(int size) {
       if (0 <= size < 10000) {
	   out_array = new out_class[size];
	 }
       else{ fprintf(outfile,"\nWARNING: bad number of out of plane angles\n");
	 }
     }

   ~out_set() {
     delete[] out_array;
     }

    void set_num(int i) { num = i;}
    int  get_num(void) { return num;}
    void set_id(int index, int new_id) { out_array[index].set_id(new_id);}
    int  get_id(int index) { return out_array[index].get_id();}
    void set_A(int index, int new_A) { out_array[index].set_A(new_A);}
    int  get_A(int index) {return out_array[index].get_A();}
    void set_B(int index, int new_B) { out_array[index].set_B(new_B);}
    int  get_B(int index) { return out_array[index].get_B();}
    void set_C(int index, int new_C) { out_array[index].set_C(new_C);}
    int  get_C(int index) { return out_array[index].get_C();}
    void set_D(int index, int new_D) { out_array[index].set_D(new_D);}
    int  get_D(int index) { return out_array[index].get_D();}
    void set_val(int index, double new_val) { out_array[index].set_value(new_val);}
    double  get_val(int index) { return out_array[index].get_value();}

    void set_s_A(int index, double s_A0, double s_A1, double s_A2) {
                 out_array[index].set_s_A(s_A0,s_A1,s_A2); }
    double get_s_A(int index, int i) { return out_array[index].get_s_A(i); }

    void set_s_B(int index, double s_B0, double s_B1, double s_B2) {
                 out_array[index].set_s_B(s_B0,s_B1,s_B2); }
    double get_s_B(int index, int i) { return out_array[index].get_s_B(i); }

    void set_s_C(int index, double s_C0, double s_C1, double s_C2) {
                 out_array[index].set_s_C(s_C0,s_C1,s_C2); }
    double get_s_C(int index, int i) { return out_array[index].get_s_C(i); }

    void set_s_D(int index, double s_D0, double s_D1, double s_D2) {
                 out_array[index].set_s_D(s_D0,s_D1,s_D2); }
    double get_s_D(int index, int i) { return out_array[index].get_s_D(i); }
    void print(FILE *fp_out, int print_flag) {
      int i;
      if (num > 0) {
        if (print_flag == 0) fprintf(fp_out,"  out = (\n");
        else fprintf (fp_out,"Out-of-planes\n");
        for (i=0; i < num; ++i)
           out_array[i].print(fp_out, print_flag);
        if (print_flag == 0) fprintf(fp_out,"  )\n");
      }
      return;
    }
    void compute(int num_atoms, double *geom) {
      int i,j,A,B,C,D;
      double rBA, rBC, rBD, phi_CBD = 0.0, dotprod = 0.0, angle = 0.0;
      double eBA[3], eBC[3], eBD[3], tmp[3];
      double *geom_ang;
    
      geom_ang = init_array(num_atoms*3);
      for (i=0;i<num_atoms*3;++i)
        geom_ang[i] = geom[i] * _bohr2angstroms;
    
      for (i=0;i<num;++i) {
        A = get_A(i);
        B = get_B(i);
        C = get_C(i);
        D = get_D(i);
    
        for (j=0;j<3;++j) {
          eBA[j] = geom_ang[3*A+j] - geom_ang[3*B+j];
          eBC[j] = geom_ang[3*C+j] - geom_ang[3*B+j];
          eBD[j] = geom_ang[3*D+j] - geom_ang[3*B+j];
        }
    
        rBA = sqrt( SQR(eBA[0]) + SQR(eBA[1]) + SQR(eBA[2]) );
        rBC = sqrt( SQR(eBC[0]) + SQR(eBC[1]) + SQR(eBC[2]) );
        rBD = sqrt( SQR(eBD[0]) + SQR(eBD[1]) + SQR(eBD[2]) );
    
        scalar_div(rBA,eBA);
        scalar_div(rBC,eBC);
        scalar_div(rBD,eBD);
    
        dot_arr(eBC,eBD,3,&phi_CBD);
    
        if (phi_CBD > 1.0) phi_CBD = 0.0;
        else if (phi_CBD < -1.0) phi_CBD = _pi ;
        else phi_CBD = acos(phi_CBD) ;
    
        cross_product(eBC,eBD,tmp);
    
        dot_arr(tmp,eBA,3,&dotprod);
    
        if (sin(phi_CBD) > optinfo.sin_phi_denominator_tol) dotprod = dotprod / sin(phi_CBD) ;
        else dotprod = 0.0 ;
    
        if (dotprod > 1.0) angle = _pi / 2.0;
        else if (dotprod < -1.0) angle = -1.0 * _pi / 2.0000;
        else angle = asin(dotprod) ;
    
        set_val(i,angle*180.0/_pi);
      }
      return;
    }
    void compute_s(int num_atoms, double *geom) {
      int i,j,A,B,C,D;
      double rBA, rBC, rBD, phi_CBD = 0.0, val = 0.0;
      double eBA[3], eBC[3], eBD[3], tmp[3], tmp2[3], tmp3[3], temp;
      double *geom_ang;
    
      geom_ang = init_array(num_atoms*3);
      for (i=0;i<num_atoms*3;++i)
        geom_ang[i] = geom[i] * _bohr2angstroms;
    
      for (i=0;i<num;++i) {
        A = get_A(i);
        B = get_B(i);
        C = get_C(i);
        D = get_D(i);
        val = get_val(i)*_pi/180.0;
    
    //fprintf(outfile,"val: %15.10lf\n",val);
    
        for (j=0;j<3;++j) {
          eBA[j] = geom_ang[3*A+j] - geom_ang[3*B+j];
          eBC[j] = geom_ang[3*C+j] - geom_ang[3*B+j];
          eBD[j] = geom_ang[3*D+j] - geom_ang[3*B+j];
        }
    
        rBA = sqrt( SQR(eBA[0]) + SQR(eBA[1]) + SQR(eBA[2]) );
        rBC = sqrt( SQR(eBC[0]) + SQR(eBC[1]) + SQR(eBC[2]) );
        rBD = sqrt( SQR(eBD[0]) + SQR(eBD[1]) + SQR(eBD[2]) );
    
    //fprintf(outfile,"rBA: %15.10lf, rBC: %15.10lf, rBD:%15.10lf\n",rBA,rBC,rBD);
    
        scalar_div(rBA,eBA);
        scalar_div(rBC,eBC);
        scalar_div(rBD,eBD);
    
    //fprintf(outfile,"eBA: %15.10lf %15.10lf %15.10lf\n",eBA[0],eBA[1],eBA[2]);
    //fprintf(outfile,"eBC: %15.10lf %15.10lf %15.10lf\n",eBC[0],eBC[1],eBC[2]);
    //fprintf(outfile,"eBD: %15.10lf %15.10lf %15.10lf\n",eBD[0],eBD[1],eBD[2]);
    
        dot_arr(eBC,eBD,3,&phi_CBD);
    
        if (phi_CBD > 1.0) phi_CBD = 0.0;
        else if (phi_CBD < -1.0) phi_CBD = _pi;
        else phi_CBD = acos(phi_CBD);
    
    //fprintf(outfile,"phi_CBD: %15.10lf\n",phi_CBD);
    
        cross_product(eBC,eBD,tmp);
        scalar_div(cos(val)*sin(phi_CBD),tmp);
        for (j=0;j<3;++j) 
           tmp2[j] = tan(val) * eBA[j];
        for (j=0;j<3;++j) 
           tmp3[j] = (tmp[j] - tmp2[j])/rBA;
        set_s_A(i,tmp3[0],tmp3[1],tmp3[2]);
    
        cross_product(eBD,eBA,tmp);
        scalar_div(cos(val)*sin(phi_CBD),tmp);
        for (j=0;j<3;++j)
          tmp2[j] = cos(phi_CBD) * eBD[j];
        for (j=0;j<3;++j)
          tmp3[j] = eBC[j] - tmp2[j];
        scalar_mult(tan(val)/SQR(sin(phi_CBD)),tmp3,3);
        for (j=0;j<3;++j)
           tmp2[j] = (tmp[j] - tmp3[j])/rBC;
        set_s_C(i,tmp2[0],tmp2[1],tmp2[2]);
    
        cross_product(eBA,eBC,tmp);
        scalar_div(cos(val)*sin(phi_CBD),tmp);
        for (j=0;j<3;++j)
          tmp2[j] = cos(phi_CBD) * eBC[j];
        for (j=0;j<3;++j)
          tmp3[j] = eBD[j] - tmp2[j];
        scalar_mult(tan(val)/SQR(sin(phi_CBD)),tmp3,3);
        for (j=0;j<3;++j)
           tmp2[j] = (tmp[j] - tmp3[j])/rBD;
        set_s_D(i,tmp2[0],tmp2[1],tmp2[2]);
    
        tmp[0] = tmp[1] = tmp[2] = 0.0;
        for (j=0;j<3;++j) {
          tmp[j]  = (-1.0) * get_s_A(i,j);
          tmp[j] -= get_s_C(i,j);
          tmp[j] -= get_s_D(i,j);
        }
        set_s_B(i,tmp[0],tmp[1],tmp[2]);
    
      }
      return;
    }
    void print_s() {
      int i;
      for (i=0;i<num;++i) {
        fprintf(outfile,"S vector for out %d %d %d %d:atom A\n",
                          get_A(i),get_B(i),get_C(i),get_D(i) );
        fprintf(outfile,"(%16.10f,%16.10f,%16.10f)\n",
                          get_s_A(i,0), get_s_A(i,1), get_s_A(i,2) );
        fprintf(outfile,"S vector for out %d %d %d %d:atom B\n",
                          get_A(i),get_B(i),get_C(i),get_D(i) );
            fprintf(outfile,"(%16.10f,%16.10f,%16.10f)\n",
                          get_s_B(i,0), get_s_B(i,1), get_s_B(i,2) );
        fprintf(outfile,"S vector for out %d %d %d %d: atom C\n",
                              get_A(i),get_B(i),get_C(i),get_D(i) );
        fprintf(outfile,"(%16.10f,%16.10f,%16.10f)\n",
                          get_s_C(i,0), get_s_C(i,1), get_s_C(i,2) );
        fprintf(outfile,"S vector for out %d %d %d %d: atom D\n",
                          get_A(i),get_B(i),get_C(i),get_D(i) );
        fprintf(outfile,"(%16.10f,%16.10f,%16.10f)\n",
                          get_s_D(i,0), get_s_D(i,1), get_s_D(i,2) );
      }
      return;
    }
    int get_id_from_atoms(int a, int b, int c, int d, int *sign) {
       *sign = 1;
       int i,A,B,C,D;
       for (i=0;i<num;++i) {
         A = get_A(i);
         B = get_B(i);
         C = get_C(i);
         D = get_D(i);
         if ( (a == A) && (b == B) && (c == C) && (d == D)) break;
         if ( (a == A) && (b == B) && (c == D) && (d == C)) {
           *sign = -1;
fprintf(outfile,"ordering change of out of plane %d\n",i);
           break;
         }
       }
       return get_id(i);
    }
};

