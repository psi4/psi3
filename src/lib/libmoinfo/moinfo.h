#ifndef _psi_src_lib_libmoinfo_moinfo_h_
#define _psi_src_lib_libmoinfo_moinfo_h_

#include <bitset>
#include <string>
#include <vector>
#include <utility>

#include "moinfo_base.h"


#define size_det 2048

#define MRCC_ON_DISK 100
#define MRCC_SO_INTS 101

enum spin                 {alpha, beta};
enum scftype              {rhf,uhf,rohf,tcscf};

class MOInfo : public MOInfoBase
{
  typedef std::vector<std::string>            strvec;
  typedef std::vector<int>                    intvec;
  typedef std::vector<std::pair<int,int> >    intpairvec;
public:
  /*********************************************************
    SlaterDeterminant Class
    1) Purpose
      This class is used to store all the information that
      belongs to a Slater Determinant
    2) Use
    3) Details
      The MOs that describe the reference are stored in the
      arrays aocc,bocc,avir and bvir. These refer to the MOs
      in Pitzer order with the frozen occupied and virtual
      already excluded. Therefore this code assumes that the
      integral transformation code has already eliminated
      the frozen integrals and relabeled them.

      type stores the type of reference
      type = 0 -> closed-shell determinant
      type = 2 -> open-shell   determinant

      number stores the ID of this reference

    4) Uses
      STL vector
  *********************************************************/
  class SlaterDeterminant{
  public:
    typedef std::bitset<size_det> bitdet;
    SlaterDeterminant();
    ~SlaterDeterminant();
    void        print();
    void        print_occ();
    void        set(int n)                               {bits.set(n);}
    bool        test(int n)                        const {return(bits.test(n));}
    bool        is_closed_shell();
    bool        is_spin_flipped(SlaterDeterminant& det);
    bitdet&     get_bits()                               {return(bits);}
    void        get_internal_excitations(SlaterDeterminant& det, double& sign,
                                                  std::vector<std::pair<int,int> >& alpha_operators,
                                                  std::vector<std::pair<int,int> >& beta_operators);
    char        get_occupation_symbol(int i);
    intvec      get_aocc();
    intvec      get_bocc();
    intvec      get_avir();
    intvec      get_bvir();

  private:
    double      annihilate(bitdet& bits_det,int so);
    double      create(bitdet& bits_det,int so);
    int         range;
    bitdet      bits;
    std::string type;
  };



public:
  friend class SlaterDeterminant;
  MOInfo();
  ~MOInfo();

  // DGEMM timing
  void        set_dgemm_timing(double value)           {dgemm_timing=value;}
  void        add_dgemm_timing(double value)           {dgemm_timing+=value;}
  double      get_dgemm_timing()                 const {return(dgemm_timing);}

  // Convergency Options
  double      get_no_damp_convergence()          const {return(no_damp_convergence);}




  int         get_mo_sym(int i)                  const {return(mo_irr[i]);}

  int         get_root()                         const {return(root);}
  

  int         get_nmo()                          const {return(nmo);}
  int         get_norbs()                        const {return(norbs);}
  int         get_nactive_ael()                  const {return(nactive_ael);}
  int         get_nactive_bel()                  const {return(nactive_bel);}
  int         get_nael()                         const {return(nael);}
  int         get_nbel()                         const {return(nbel);}

  int         get_nfocc()                        const {return(nfocc);}
  int         get_navir()                        const {return(navir);}
  int         get_nfvir()                        const {return(nfvir);}
  int         get_nocc()                         const {return(nocc);}
  int         get_nvir()                         const {return(nvir);}


  int*        get_orbspi()                       const {return(orbspi);}
  int*        get_focc()                         const {return(focc);}

  int*        get_sopi()                         const {return(sopi);}
  int*        get_docc()                         const {return(docc);}
  int*        get_actv()                         const {return(actv);}
  
  int*        get_avir()                         const {return(avir);}
  int*        get_fvir()                         const {return(fvir);}
  int*        get_occ()                          const {return(occ);}
  int*        get_vir()                          const {return(vir);}

  int         get_sopi(int i)                    const {return(sopi[i]);}
  int         get_orbspi(int i)                  const {return(orbspi[i]);}
  int         get_focc(int i)                    const {return(focc[i]);}
  int         get_docc(int i)                    const {return(docc[i]);}
  int         get_actv(int i)                    const {return(actv[i]);}
  int         get_avir(int i)                    const {return(avir[i]);}
  int         get_fvir(int i)                    const {return(fvir[i]);}

  int*        get_clsdpi()                       const {return(clsdpi);}
  int*        get_openpi()                       const {return(openpi);}

  // Mapping functions
  int         get_nonfrozen_to_all(int i)        const {return(nonfrozen_to_all[i]);}
  int         get_all_to_nonfrozen(int i)        const {return(all_to_nonfrozen[i]);}
  int*        get_first_so_pitzer()              const {return(first_so_pitzer);}
  int*        get_last_so_pitzer()               const {return(last_so_pitzer);}
  int*        get_first_orbs_pitzer()            const {return(first_orbs_pitzer);}
  int*        get_last_orbs_pitzer()             const {return(last_orbs_pitzer);}
  int         get_first_orbs_pitzer(int i)       const {return(first_orbs_pitzer[i]);}
  int         get_last_orbs_pitzer(int i)        const {return(last_orbs_pitzer[i]);}
  int*        get_first_occupied_pitzer(int i)   const {return(first_occupied_pitzer[i]);}
  int*        get_first_active_pitzer(int i)     const {return(first_active_pitzer[i]);}
  int*        get_first_virtual_pitzer(int i)    const {return(first_virtual_pitzer[i]);}
  int*        get_last_occupied_pitzer(int i)    const {return(last_occupied_pitzer[i]);}
  int*        get_last_active_pitzer(int i)      const {return(last_active_pitzer[i]);}
  int*        get_last_virtual_pitzer(int i)     const {return(last_virtual_pitzer[i]);}


  int*        get_so_to_pitzer()                 const {return(so_to_pitzer);}
  int*        get_orbs_to_pitzer()               const {return(orbs_to_pitzer);}
  int*        get_docc_to_pitzer()               const {return(docc_to_pitzer);}
  int*        get_act_to_pitzer()                const {return(act_to_pitzer);}
  int*        get_ext_to_pitzer()                const {return(ext_to_pitzer);}
  int*        get_occ_to_pitzer()                const {return(occ_to_pitzer);}
  int*        get_vir_to_pitzer()                const {return(vir_to_pitzer);}
  int*        get_all_to_pitzer()                const {return(all_to_pitzer);}
  int*        get_qt_to_pitzer()                 const {return(all_to_pitzer);}
  int*        get_actv_to_occ()                  const {return(act_to_occ);}
  int*        get_actv_to_vir()                  const {return(act_to_vir);}
  int*        get_occ_to_actv()                  const {return(occ_to_act);}
  int*        get_vir_to_actv()                  const {return(vir_to_act);}
  bool*       get_is_act_in_occ()                const {return(is_act_in_occ);}
  bool*       get_is_act_in_vir()                const {return(is_act_in_vir);}
//   vector<int> get_mapping(char* ,"pitzer")        const {return(all_to_pitzer);}



  int         get_all_to_occ(int i)              const {return(all_to_occ[i]);}
  int         get_all_to_vir(int i)              const {return(all_to_vir[i]);}
  int         get_all_to_pitzer(int i)           const {return(all_to_pitzer[i]);}

  double*     get_evals(spin s)                  const {return(evals[s]);}

  double      get_scf_energy()                   const {return(scf_energy);}
  double      get_fzcore_energy()                const {return(fzcore_energy);}
  void        set_fzcore_energy(double efzc)           {fzcore_energy=efzc;}

  // Model space functions
  void        setup_model_space();
  int         get_nrefs()                              {return(all_refs.size());};
  int         get_nunique()                            {return(unique_refs.size());};
  int         get_ref_number(std::string str,int n);
  int         get_ref_size(std::string str);
  strvec      get_matrix_names(std::string str);
  intvec      get_aocc(std::string str,int i);
  intvec      get_bocc(std::string str,int i);
  intvec      get_avir(std::string str,int i);
  intvec      get_bvir(std::string str,int i);
  intpairvec  get_alpha_internal_excitation(int i,int j);
  intpairvec  get_beta_internal_excitation(int i,int j);
  double      get_sign_internal_excitation(int i,int j);

private:
  void        tuning();
  void        read_info();
  void        read_mo_spaces();
  void        compute_mo_mappings();
  void        print_info();
  void        print_mo();
  void        free_memory_info();
  void        free_memory_mo_spaces();

  // Model space functions
  void        print_model_space();
  void        build_model_space();
  void        make_internal_excitations();

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // MOInfo variables
  int         root;


  double      scf_energy;
  double      fzcore_energy;

  /////////////////////////////////////////////////////////////////////////////////////////////////



  double      dgemm_timing;

  // In-core/Out-of-core
//   double      block_size;
  double      no_damp_convergence;

  int         nel;

  int         reference;
  // Total number of orbitals in each space

  int         norbs;      // Psi nmo

  int         nfocc;

  int         navir;
  int         nfvir;
  int         nactv_docc;
  int         nocc;                                  // Generalized occupied (docc + actv)
  int         nvir;                                  // Generalized virtual (actv + ext)
  // Orbitals arrays, should probably turn this into vector<int> soon.
  int*        focc;


  int*        avir;
  int*        fvir;
  int*        occ;
  int*        vir;
  int*        actv_docc;
  double*     evals[2];                              // Fock matrix diagonal elements
  int*        clsdpi;
  int*        openpi;

  int*        orbspi;

  // Mapping arrays 

  int*        nonfrozen_to_all;
  int*        all_to_nonfrozen;
  int*        so_to_pitzer;
  int*        orbs_to_pitzer;
  int*        docc_to_pitzer;
  int*        act_to_pitzer;
  int*        ext_to_pitzer;
  int*        occ_to_pitzer;
  int*        vir_to_pitzer;
  int*        all_to_pitzer;
  int*        pitzer_to_occ_act;
  int*        pitzer_to_act_vir;
  int*        occ_to_vir;
  int*        all_to_occ;
  int*        all_to_vir;
  int*        act_to_occ;
  int*        act_to_vir;
  int*        occ_to_act;
  int*        vir_to_act;
  bool*       is_act_in_occ;
  bool*       is_act_in_vir;
  // First-last arrays
  int*        first_orbs_pitzer;
  int*        last_orbs_pitzer;
  int*        first_so_pitzer;
  int*        last_so_pitzer;
  int*        first_occupied_pitzer[2];
  int*        first_virtual_pitzer[2];
  int*        first_active_pitzer[2];
  int*        last_occupied_pitzer[2];
  int*        last_active_pitzer[2];
  int*        last_virtual_pitzer[2];
  // Symmetry
  int*        mo_irr;

  // Model space
  std::vector<SlaterDeterminant> references;
  std::vector<std::vector<std::vector<std::pair<int,int> > > > alpha_internal_excitations;
  std::vector<std::vector<std::vector<std::pair<int,int> > > > beta_internal_excitations;
  std::vector<std::vector<double> >                            sign_internal_excitations;
  std::vector<int> all_refs;
  std::vector<int> unique_refs;
  std::vector<int> closed_shell_refs;
  std::vector<int> unique_open_shell_refs;
};

extern MOInfo  *moinfo;

#endif // _psi_src_lib_libmoinfo_moinfo_h_
