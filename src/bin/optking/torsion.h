/*! \file 
    \ingroup (OPTKING)
    \brief Enter brief description of file here 
*/
class torsion_class {
    int id;
    int A;
    int B;
    int C;
    int D;
    double value;
    double *s_A; // The s vector for atom A
    double *s_B;
    double *s_C;
    double *s_D;
    int near_lin;
    // +1 if approaching 180
    //  0 if OK
    // -1 if approaching -180
  public:
    torsion_class(){
      s_A = new double [3];
      s_B = new double [3];
      s_C = new double [3];
      s_D = new double [3];
      near_lin = 0;
    }
    ~torsion_class(){
      // fprintf(stdout,"destructing torsion class\n");
      delete [] s_A;
      delete [] s_B;
      delete [] s_C;
      delete [] s_D;
    }
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
    void set_near_lin(int new_near_lin) { near_lin = new_near_lin;}
    double get_value(void)  { return value;}
    int get_near_lin(void)  { return near_lin;}
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

class torsion_set {

  int num;
   torsion_class *tors_array;

  public:

   torsion_set(int size){
     if(0 <= size < 10000)
	   tors_array = new torsion_class[size];
     else 
       fprintf(outfile,"\nWARNING: bad number of torsions\n");
   }

   torsion_set(void){ } /* don't allocate memory yet */

   void allocate(int size) {
     if(0 <= size < 10000)
       tors_array = new torsion_class[size];
     else
       fprintf(outfile,"\nWARNING: bad number of torsions\n");
   }

   ~torsion_set() {
       // fprintf(stdout,"destructing torsion_set\n");
       delete [] tors_array;
   }

   void print(FILE *fp_out, int print_flag) {
      int i;
      if (num > 0) {
        if (print_flag == 0) fprintf(fp_out,"  tors = (\n");
        else fprintf(fp_out,"Torsions\n");
        for (i=0; i < num; ++i)
           tors_array[i].print(fp_out, print_flag);
        if (print_flag == 0) fprintf(fp_out,"  )\n");
      }
      return;
    }

   void print_s() {  int i;
      for (i=0;i<num;++i) {
        fprintf(outfile,"S vector for tors %d %d %d %d:atom A\n",
                          get_A(i),get_B(i),get_C(i),get_D(i) );
        fprintf(outfile,"(%16.10f,%16.10f,%16.10f)\n",
                          get_s_A(i,0), get_s_A(i,1), get_s_A(i,2) );
        fprintf(outfile,"S vector for tors %d %d %d %d:atom B\n",
                          get_A(i),get_B(i),get_C(i),get_D(i) );
        fprintf(outfile,"(%16.10f,%16.10f,%16.10f)\n",
                          get_s_B(i,0), get_s_B(i,1), get_s_B(i,2) );
        fprintf(outfile,"S vector for tors %d %d %d %d: atom C\n",
                          get_A(i),get_B(i),get_C(i),get_D(i) );
        fprintf(outfile,"(%16.10f,%16.10f,%16.10f)\n",
                          get_s_C(i,0), get_s_C(i,1), get_s_C(i,2) );
        fprintf(outfile,"S vector for tors %d %d %d %d: atom D\n",
                          get_A(i),get_B(i),get_C(i),get_D(i) );
        fprintf(outfile,"(%16.10f,%16.10f,%16.10f)\n",
                          get_s_D(i,0), get_s_D(i,1), get_s_D(i,2) );
      }
      return;
    }

    void set_num(int i) { num = i;}
    int  get_num(void) { return num;}
    void set_id(int index, int new_id) { tors_array[index].set_id(new_id);}
    int  get_id(int index) { return tors_array[index].get_id();}
    void set_A(int index, int new_A) { tors_array[index].set_A(new_A);}
    int  get_A(int index) {return tors_array[index].get_A();}
    void set_B(int index, int new_B) { tors_array[index].set_B(new_B);}
    int  get_B(int index) { return tors_array[index].get_B();}
    void set_C(int index, int new_C) { tors_array[index].set_C(new_C);}
    int  get_C(int index) { return tors_array[index].get_C();}
    void set_D(int index, int new_D) { tors_array[index].set_D(new_D);}
    int  get_D(int index) { return tors_array[index].get_D();}
    void set_val(int index, double new_val) { tors_array[index].set_value(new_val);}
    double  get_val(int index) { return tors_array[index].get_value();}
    void set_near_lin(int index, int new_near_lin)
        { tors_array[index].set_near_lin(new_near_lin);}
    int get_near_lin(int index)
        { return tors_array[index].get_near_lin();}
    void set_s_A(int index, double s_A0, double s_A1, double s_A2) {
                 tors_array[index].set_s_A(s_A0,s_A1,s_A2); }
    double get_s_A(int index, int i) { return tors_array[index].get_s_A(i); }

    void set_s_B(int index, double s_B0, double s_B1, double s_B2) {
                 tors_array[index].set_s_B(s_B0,s_B1,s_B2); }
    double get_s_B(int index, int i) { return tors_array[index].get_s_B(i); }

    void set_s_C(int index, double s_C0, double s_C1, double s_C2) {
                 tors_array[index].set_s_C(s_C0,s_C1,s_C2); }
    double get_s_C(int index, int i) { return tors_array[index].get_s_C(i); }

    void set_s_D(int index, double s_D0, double s_D1, double s_D2) {
                 tors_array[index].set_s_D(s_D0,s_D1,s_D2); }
    double get_s_D(int index, int i) { return tors_array[index].get_s_D(i); }

    void compute(int natom, double *geom) {
      int i,j,k,A,B,C,D,sign;
      double rAB,rBC,rCD,phi_123,phi_234,val = 0.0;
      double eAB[3], eBC[3], eCD[3], tmp[3], tmp2[3], tmp3[3];
      double *geom_ang, dotprod, angle;
    
      geom_ang = new double [3*natom];
      for (i=0;i<natom*3;++i)
        geom_ang[i] = geom[i] * _bohr2angstroms;
    
      for (i=0;i<num;++i) {
        A = get_A(i);
        B = get_B(i);
        C = get_C(i);
        D = get_D(i);
    
        for (j=0;j<3;++j) {
          eAB[j] = geom_ang[3*B+j] - geom_ang[3*A+j];
          eBC[j] = geom_ang[3*C+j] - geom_ang[3*B+j];
          eCD[j] = geom_ang[3*D+j] - geom_ang[3*C+j];
        }
    
        rAB = sqrt( SQR(eAB[0]) + SQR(eAB[1]) + SQR(eAB[2]) );
        rBC = sqrt( SQR(eBC[0]) + SQR(eBC[1]) + SQR(eBC[2]) );
        rCD = sqrt( SQR(eCD[0]) + SQR(eCD[1]) + SQR(eCD[2]) );
    
        scalar_div(rAB,eAB);
        scalar_div(rBC,eBC);
        scalar_div(rCD,eCD);
    
    /*
      fprintf(outfile,"eAB\n");
      fprintf(outfile,"%20.16lf %20.16lf %20.16lf \n",eAB[0],eAB[1],eAB[2]);
      fprintf(outfile,"eBC\n");
      fprintf(outfile,"%20.16lf %20.16lf %20.16lf \n",eBC[0],eBC[1],eBC[2]);
      fprintf(outfile,"eCD\n");
      fprintf(outfile,"%20.16lf %20.16lf %20.16lf \n",eCD[0],eCD[1],eCD[2]);
    */
    
        phi_123 = 0.0;
        phi_234 = 0.0;
    
        for (j=0;j<3;++j) {
           phi_123 += (-1.0*eAB[j]) * eBC[j];
           phi_234 += (-1.0*eBC[j]) * eCD[j];
        }
    
        if (phi_123 > 1.0)
          phi_123 = 0.0;
        else if (phi_123 < -1.0)
          phi_123 = _pi;
        else phi_123 = acos(phi_123);
    
        if (phi_234 > 1.0)
          phi_234 = 0.0;
        else if (phi_234 < -1.0)
          phi_234 = _pi;
        else phi_234 = acos(phi_234);
    
        cross_product(eAB,eBC,tmp);
        cross_product(eBC,eCD,tmp2);
        dot_arr(tmp,tmp2,3,&dotprod);

        if ((sin(phi_123) > optinfo.sin_phi_denominator_tol) &&
            (sin(phi_234) > optinfo.sin_phi_denominator_tol)) {
           dotprod /= sin(phi_123);
           dotprod /= sin(phi_234);
        }
        else dotprod = 2.0 ;
    
        if (dotprod > optinfo.cos_tors_near_1_tol) angle = 0.0000 ;
        else if (dotprod < optinfo.cos_tors_near_neg1_tol) angle = _pi ;
        else angle = acos(dotprod) ;

        // determine sign of torsions
        //cross_product(tmp,tmp2,tmp3);
        //dot_arr(tmp3,eBC,3,&dotprod);
        cross_product(eBC,eCD,tmp);
        dot_arr(eAB,tmp,3,&dotprod);
        if (dotprod < 0) sign = -1; else sign = 1;

        // extend domain of torsions so delta(values) can be calculated
        angle = sign * angle * 180.0 / _pi;
        if ((get_near_lin(i) == -1) && (angle > 160.0)) {
// fprintf(outfile,"get_near_lin(%d)=%d angle %15.10lf angle %15.10lf\n",
    // i, get_near_lin(i), angle, -180.0 - (180.0 - angle) );
          angle = -180.0 - (180.0 - angle);
        }
        else if ((get_near_lin(i) == +1) && (angle < -160.0)) {
// fprintf(outfile,"get_near_lin(%d)=%d angle %15.10lf angle %15.10lf\n",
    // i, get_near_lin(i), angle, +180.0 + (180.0 + angle) );
          angle = +180.0 + (180.0 + angle);
        }

        set_val(i,angle);
      }
      delete [] geom_ang;
      return;
    }
    void fix_near_lin(void) {
      int i, lin;
      for (i=0;i<num;++i) {
        if ( get_val(i) > 160.0) {
          lin = +1;
        }
        else if ( get_val(i) < -160.0) {
          lin = -1;
        }
        else
          lin = 0;
        set_near_lin(i,lin);
      }
      return;
    }
    void compute_s(int natom, double *geom) {
      int i,j,A,B,C,D;
      double rAB,rBC,rCD;
      double eAB[3], eBC[3], eCD[3], tmp[3], tmp2[3];
      double phiABC, phiBCD;
      double *geom_ang;
    
      geom_ang = new double [3*natom];
      for (i=0;i<natom*3;++i)
        geom_ang[i] = geom[i] * _bohr2angstroms;
    
      for (i=0;i<num;++i) {
        A = get_A(i);
        B = get_B(i);
        C = get_C(i);
        D = get_D(i);
    
        for (j=0;j<3;++j) {
          eAB[j] = geom_ang[3*B+j] - geom_ang[3*A+j];
          eBC[j] = geom_ang[3*C+j] - geom_ang[3*B+j];
          eCD[j] = geom_ang[3*D+j] - geom_ang[3*C+j];
        }
    
        rAB = sqrt( SQR(eAB[0]) + SQR(eAB[1]) + SQR(eAB[2]) );
        rBC = sqrt( SQR(eBC[0]) + SQR(eBC[1]) + SQR(eBC[2]) );
        rCD = sqrt( SQR(eCD[0]) + SQR(eCD[1]) + SQR(eCD[2]) );
    
        scalar_div(rAB,eAB);
        scalar_div(rBC,eBC);
        scalar_div(rCD,eCD);
    
        phiABC = 0.0;
        phiBCD = 0.0;
    
        for (j=0;j<3;++j) {
           phiABC += (-1.0 * eAB[j]) * eBC[j];
           phiBCD += (-1.0 * eBC[j]) * eCD[j];
        }
    
        phiABC = acos(phiABC);
        phiBCD = acos(phiBCD);
    
        cross_product(eAB,eBC,tmp);
        scalar_div(-1.0 * rAB * SQR(sin(phiABC)),tmp);
        set_s_A(i,tmp[0],tmp[1],tmp[2]);
    
        cross_product(eAB,eBC,tmp);
        scalar_mult((rBC-rAB*cos(phiABC))/(rBC*rAB*SQR(sin(phiABC))),tmp,3);
        cross_product(eCD,eBC,tmp2);
        scalar_mult(cos(phiBCD)/(rBC*SQR(sin(phiBCD))),tmp2,3);
        set_s_B(i,tmp[0]+tmp2[0],tmp[1]+tmp2[1],tmp[2]+tmp2[2]);
    
        cross_product(eCD,eBC,tmp);
        scalar_mult((rBC-rCD*cos(phiBCD))/(rBC*rCD*SQR(sin(phiBCD))),tmp,3);
        cross_product(eAB,eBC,tmp2);
        scalar_mult(cos(phiABC)/(rBC*SQR(sin(phiABC))),tmp2,3);
        set_s_C(i,tmp[0]+tmp2[0],tmp[1]+tmp2[1],tmp[2]+tmp2[2]);
    
        cross_product(eCD,eBC,tmp);
        scalar_div(-1.0*rCD*SQR(sin(phiBCD)),tmp);
        set_s_D(i,tmp[0],tmp[1],tmp[2]);
      }
      delete [] geom_ang;
      return;
    }
    int get_id_from_atoms(int a, int b, int c, int d) {
       int i, foundit = 0;
//   fprintf(outfile,"tors.get_id_from_atoms(%d,%d,%d,%d)\n",a,b,c,d);
       for (i=0;i<num;++i) {
         if ( (a == get_A(i)) && (b == get_B(i))
           && (c == get_C(i)) && (d == get_D(i)))
             return get_id(i);
       }

       for (i=0;i<num;++i) {
         if ( (a == get_A(i)) && (c == get_B(i))
           && (b == get_C(i)) && (d == get_D(i)))
             return get_id(i);
       }
       if (optinfo.delocalize) {
         if (i == num) {
           fprintf(outfile,"Could not find simple torsion for atoms  \
               %d %d %d %d in list.\n", a+1, b+1, c+1, d+1);
           exit(2);
         }
       }
       return -1;
//   fprintf(outfile,"Returning id: %d\n",get_id(i));
     }
};

