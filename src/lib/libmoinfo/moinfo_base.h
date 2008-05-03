#ifndef _psi_src_lib_libmoinfo_moinfo_base_h_
#define _psi_src_lib_libmoinfo_moinfo_base_h_

/*! \file 
    \ingroup LIBMOINFO
    \brief   This class stores all the basic info regarding MOs
*/

#define PSI_NULL(args) args = NULL;
#define PSI_FREE(args) if(args != NULL) free(args);
#define PSI_DELETE(args) if(args != NULL) delete args;
#define PSI_DELETE_ARRAY(args) if(args != NULL) delete[] args;
#define IOFF 500000

#include <string>

#include <libutil/libutil.h>

class MOInfoBase{
public:
  MOInfoBase();
  ~MOInfoBase();
  
  double      get_nuclear_energy()               const {return(nuclear_energy);}
  
  char**      get_irr_labs()                     const {return(irr_labs);}
  char*       get_irr_labs(int i)                const {return(irr_labs[i]);}
  
  int         get_nirreps()                      const {return(nirreps);}
  int         get_nso()                          const {return(nso);}
  
  int*        get_ioff()                         const {return(ioff);}
  int*        get_sopi()                         const {return(sopi);}
  int*        get_docc()                         const {return(docc);}
  int*        get_actv()                         const {return(actv);}
  
  int         get_ndocc()                        const {return(ndocc);}
  int         get_nactv()                        const {return(nactv);}
  
  double**    get_scf_mos()                      const {return(scf);}
  double**    get_scf_mos(int i)                 const {return(scf_irrep[i]);}
  double      get_scf_mos(int i,int j)           const {if((i<nmo)&&(j<nso)) return(scf[i][j]); else print_error("get_scf_mos out of range",__FILE__,__LINE__); return(0.0);}
  void        write_chkpt_mos();
protected:
  void        read_chkpt_data();
  void        compute_number_of_electrons();
  void        correlate(char *ptgrp, int irrep, int& nirreps_old, int& nirreps_new,int*& correlation);
  void        read_mo_space(int nirreps_ref,int& n, int* mo, const char* label);
  void        print_mo_space(int& nmo, int* mo, const char* label);
  
  void        startup();
  void        cleanup();
  void        compute_ioff();

  int         nirreps;
  int         wfn_sym;
  int         charge;
  int         multiplicity;
  
  int         nso;              // PSI nso (number of symmetry-adapted atomic orbitals)
  int         nmo;              // Psi nmo (number of molecular orbitals, including frozen core and frozen virtual)
  int         ndocc;
  int         nactv;
  int         nael;
  int         nbel;
  int         nactive_ael;
  int         nactive_bel;
  
  int*        ioff;
  int*        sopi;
  int*        docc;
  int*        actv;
  
  double      nuclear_energy;
  
  double**    scf;                                   // MO coefficients
  double***   scf_irrep;                             // MO coefficients
  
  char**      irr_labs;
};

#endif // _psi_src_lib_libmoinfo_moinfo_base_h_

/*

    MOInfoBase();
  MOInfoBase(std::string moinfo_prefix);
  ~MOInfoBase();
  MOInfoBase();
  MOInfoBase(std::string moinfo_prefix);
  ~MOInfoBase();
  
  void        print();

  // Basic PSI variables
  int         get_nirreps()                      const {return(nirreps);}
  int         get_nso()                          const {return(nso);}
  int         get_nael()                         const {return(nael);}
  int         get_nbel()                         const {return(nbel);}
  
  int*        get_sopi()                         const {return(sopi);}
  int         get_sopi(int h)                    const {return(sopi[h]);}

  int         get_ndocc()                        const {return(ndocc);}
  int         get_nactv()                        const {return(nactv);}

  int         get_docc(int h)                    const {return(docc[h]);}
  int         get_actv(int h)                    const {return(actv[h]);}

protected:
  std::string moinfo_prefix_;
  

  int         nirreps;
  int         nso;
  int         nmo;
  int         nael;
  int         nbel;
  int         nactive_ael;
  int         nactive_bel;
  int         ndocc;
  int         nactv;
  int         wfn_sym;
  
  int*        sopi;
  int*        docc;
  int*        actv;
  
  char**      irr_labs;

  void        startup();
  void        cleanup();
  void        read_chkpt_data();
  void        compute_number_of_electrons();
  void        read_mo_space(int nirreps_ref,int& n, int* mo, char* label);
  void        correlate(char *ptgrp, int irrep, int& nirreps_old, int& nirreps_new,int*& arr);*/