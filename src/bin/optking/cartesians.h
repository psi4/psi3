/*! \file
    \ingroup OPTKING
    \brief Class declaration for cartesian coordinates
*/

#ifndef _psi3_bin_optking_cartesians_h_
#define _psi3_bin_optking_cartesians_h_

namespace psi { namespace optking {

class cartesians {
    double energy;
    int natom;
    int nallatom;
    double *atomic_num;
    double *coord;
    double *grad;
    double *mass;
    double *fatomic_num;
    double *fcoord;
    double *fgrad;
    double *fmass;

  public:
    ~cartesians() {
      // fprintf(stdout,"destructing cartesian\n");
      delete [] atomic_num;
      delete [] coord;
      delete [] grad;
      delete [] mass;
      delete [] fatomic_num;
      delete [] fcoord;
      delete [] fgrad;
      delete [] fmass;
    }
    void print(int flag, FILE *fp_out, int new_geom_file, char *disp_label,
               int disp_num);

    void set_coord(double *geom) {
      int i, xyz, cnt;
      for (i=0;i<natom;++i) {
        for (xyz=0; xyz<3; ++xyz) {
          coord[3*i+xyz] = geom[3*i+xyz];
          fcoord[3*optinfo.to_dummy[i]+xyz] = geom[3*i+xyz];
        }
      }
      return;
    }

    void set_coord_2d(double **geom_2d) {
      int i, xyz, cnt=0;
      double *tmp_geom;
      tmp_geom = init_array(3*natom);
      for (i=0; i<natom; ++i)
        for (xyz=0; xyz<3; ++xyz)
          tmp_geom[cnt++] = geom_2d[i][xyz];
      set_coord(tmp_geom); 
      free(tmp_geom);
    }

    void set_fcoord(double *geom) {
      int i, cnt, xyz;
      for (i=0; i<nallatom; ++i) {
        for (xyz=0; xyz<3; ++xyz) {
          fcoord[3*i+xyz] = geom[3*i+xyz];
          coord[3*optinfo.to_nodummy[i]+xyz] = geom[3*i+xyz];
        }
      }
      return;
    }

    void set_grad(double *gradient) {
      int i;
      for (i=0;i<natom*3;++i)
        grad[i] = gradient[i];
      return;
    }
    double *get_coord() {
      int i;
      double *copy;
      copy = init_array(natom*3);
      for (i=0;i<natom*3;++i)
        copy[i] = coord[i];
      return copy;
    }
    double **get_coord_2d() {
      int i,j,cnt;
      double **copy;
      copy = block_matrix(natom,3);
      cnt=0;
      for (i=0; i<natom; ++i)
        for (j=0; j<3; ++j)
          copy[i][j] = coord[cnt++];
      return copy;
    }
    double *get_fcoord() {
      int i;
      double *copy;
      copy = init_array(nallatom*3);
      for (i=0;i<(nallatom*3);++i)
        copy[i] = fcoord[i];
      return copy;
    }
    double *get_mass() {
      int i;
      double *copy;
      copy = init_array(natom*3);
      for (i=0;i<natom*3;++i)
        copy[i] = mass[i];
      return copy;
    }
    double *get_Zvals() {
      int i;
      double *copy;
      copy = init_array(natom*3);
      for (i=0;i<natom*3;++i)
        copy[i] = atomic_num[i];
      return copy;
    }
    double *get_fmass() {
      int i;
      double *copy;
      copy = init_array(nallatom*3);
      for (i=0;i<nallatom*3;++i)
        copy[i] = fmass[i];
      return copy;
    }
    void mult(double factor) {
      int i;
      for (i=0;i<natom*3;++i) {
         coord[i] *= factor;
      }
      return;
    }
    double val(int i, int j) { return coord[3*i+j]; }
    double *get_forces();
    double *get_fforces();
    int get_natom() {return natom; }
    int get_nallatom() {return nallatom; }
    void set_natom(int new_num) {natom = new_num;}
    void set_energy(double new_energy) {energy = new_energy;}
    double get_energy() {return energy;} 
    double get_Z(int i) { return atomic_num[i]; }
    double get_fatomic_num(int i) { return fatomic_num[i]; }
    cartesians();
    double R(int i, int j) {
      int xyz;
      double tval = 0.0;
      for (xyz=0; xyz<3; ++xyz)
        tval += (coord[3*i+xyz] - coord[3*j+xyz]) * (coord[3*i+xyz] - coord[3*j+xyz]);
      return sqrt(tval);
    }
};

}} /* namespace psi::optking */

#endif
