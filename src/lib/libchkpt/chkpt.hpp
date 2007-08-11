#ifndef _psi_src_lib_libchkpt_chkpt_hpp_
#define _psi_src_lib_libchkpt_chkpt_hpp_

#include <libchkpt/config.h>

namespace psi {
	class PSIO;
	
	/**
		Chkpt is an instance of the libchkpt library. Multiple instances are 
		supported by having multiple instances of libpsio available.
		
		Each instance is configured by the constructor:
		Chkpt new_instance(&libpsio_instance, PSIO_OPEN_OLD);
	*/
	class Chkpt {
		/*! Instance of libpsio to use. */
		PSIO *psio;
		/*! Active TOC entry prefix. */
		char chkpt_prefix[CHKPT_PREFIX_LEN];
	public:
		/*! Constructor. Calls PSIO::open to open the checkpoint file.
			\param psioObject Instance of libpsio to connect through.
			\param status Either PSIO_OPEN_OLD or PSIO_OPEN_NEW.
		*/
		Chkpt(PSIO *psioObject, int status);
		/*! Destructor. Call PSIO::close to close the checkpoint file. */
		~Chkpt();

		char *rd_prefix();
		void wt_prefix(char *prefix);
		void set_prefix(char *prefix);
		void commit_prefix();
		void reset_prefix();
		char *get_prefix();

		char *build_keyword(char *key);

		int exist(char *keyword);

		char *rd_label();
		void wt_label(char *label);

		double rd_escf();
		void wt_escf(double escf);

		double rd_eref();
		void wt_eref(double eref);

		double rd_ecorr();
		void wt_ecorr(double ecorr);

		double rd_enuc();
		void wt_enuc(double enuc);

		double rd_efzc();
		void wt_efzc(double efzc);

		double rd_etot();
		void wt_etot(double etot);

		int rd_disp();
		void wt_disp(int disp);

		double rd_eccsd();
		void wt_eccsd(double eccsd);

		double rd_e_t();
		void wt_e_t(double e_t);

		double rd_emp2();
		void wt_emp2(double emp2);
		
		int *rd_am2canon_shell_order(void);
		void wt_am2canon_shell_order(int *);
		
		int* rd_atom_dummy(void);
		void wt_atom_dummy(int *);

		int *rd_atom_position(void);
		void wt_atom_position(int *);

		double **rd_cartrep();
		void wt_cartrep(double **);

		double **rd_ccvecs(void);
		void wt_ccvecs(double **);

		double** rd_cdsalc2cd();
		void wt_cdsalc2cd(const double**);

		int *rd_clsdpi(void);
		void wt_clsdpi(int *);

		double *rd_contr(void);
		void wt_contr(double *);
		double **rd_contr_full(void);

		int rd_disp_irrep(void);
		void wt_disp_irrep(int);

		double rd_e_labeled(char *);
		void wt_e_labeled(char *, double);

		double *rd_evals(void);
		double *rd_alpha_evals(void);
		double *rd_beta_evals(void);
		void wt_evals(double *);
		void wt_alpha_evals(double *);
		void wt_beta_evals(double *);

		double *rd_exps(void);
		void wt_exps(double *);

		char **rd_felement(void);
		void wt_felement(char **);

		double **rd_fgeom(void);
		void wt_fgeom(double **);

		int *rd_frzcpi(void);
		void wt_frzcpi(int *);

		int *rd_frzvpi(void);
		void wt_frzvpi(int *);

		double **rd_geom(void);
		void wt_geom(double **);

		double* rd_grad(void);
		void wt_grad(double*);		

		int **rd_ict(void);
		void wt_ict(int **);

		int rd_iopen(void);
		void wt_iopen(int);

		char **rd_irr_labs(void);
		void wt_irr_labs(char **);

		double **rd_lagr(void);
		double **rd_alpha_lagr(void);
		double **rd_beta_lagr(void);
		void wt_lagr(double **);
		void wt_alpha_lagr(double **);
		void wt_beta_lagr(double **);

		int rd_max_am(void);
		void wt_max_am(int);

		int rd_nallatom(void);
		void wt_nallatom(int);

		int rd_nao(void);
		void wt_nao(int);

		int rd_natom(void);
		void wt_natom(int);

		int rd_ncalcs(void);

		int rd_nfzc(void);
		void wt_nfzc(int);

		int rd_nfzv(void);
		void wt_nfzv(int);

		int rd_nirreps(void);
		void wt_nirreps(int);

		int rd_nmo(void);
		void wt_nmo(int);

		int rd_nprim(void);
		void wt_nprim(int);

		int rd_nshell(void);
		void wt_nshell(int);

		int rd_nso(void);
		void wt_nso(int);

		int rd_nsymhf(void);
		void wt_nsymhf(int);

		int rd_num_unique_atom(void);
		void wt_num_unique_atom(int);

		int rd_num_unique_shell(void);
		void wt_num_unique_shell(int);

		int *rd_openpi(void);
		void wt_openpi(int *);

		int *rd_orbspi(void);
		void wt_orbspi(int *);

		int rd_override_occ(void);
		void wt_override_occ(int);

		int rd_phase_check(void);
		void wt_phase_check(int);

		int rd_ref(void);
		void wt_ref(int);

		int rd_rottype(void);
		void wt_rottype(int);

		double **rd_rref(void);
		void wt_rref(double **);

		double **rd_scf(void);
		double **rd_alpha_scf(void);
		double **rd_beta_scf(void);
		void wt_scf(double **);
		void wt_alpha_scf(double **);
		void wt_beta_scf(double **);
		double **rd_local_scf(void);
		void wt_local_scf(double **);
		double **rd_scf_irrep(int irrep);
		double **rd_alpha_scf_irrep(int irrep);
		double **rd_beta_scf_irrep(int irrep);
		void wt_scf_irrep(double**, int);
		void wt_alpha_scf_irrep(double**, int);
		void wt_beta_scf_irrep(double**, int);
		double** set_mo_phases(double**, int, int);
		
		int **rd_shell_transm(void);
		void wt_shell_transm(int **);

		int *rd_sloc(void);
		void wt_sloc(int *);

		int *rd_shells_per_am(void);
		void wt_shells_per_am(int *);

		int *rd_sloc_new(void);
		void wt_sloc_new(int *);

		int *rd_snuc(void);
		void wt_snuc(int *);

		int *rd_snumg(void);
		void wt_snumg(int *);

		int *rd_sopi(void);
		void wt_sopi(int *);

		int *rd_sprim(void);
		void wt_sprim(int *);

		int *rd_statespi(void);
		void wt_statespi(int *);

		int *rd_stype(void);
		void wt_stype(int *);

		char *rd_sym_label(void);
		void wt_sym_label(char *sym_label);

		int *rd_symoper(void);
		void wt_symoper(int *);

		int *rd_ua2a(void);
		void wt_ua2a(int *);

		int *rd_us2s(void);
		void wt_us2s(int *);

		double **rd_usotao(void);
		void wt_usotao(double **);

		double **rd_usotbf(void);
		void wt_usotbf(double **);

		struct z_entry *rd_zmat(void);
		void wt_zmat(struct z_entry *);

		double *rd_zvals(void);
		void wt_zvals(double *zvals);
		
		int* rd_cdsalcpi();
		void wt_cdsalcpi(const int*);
	};
	
	extern Chkpt* _default_chkpt_lib_;
}

#endif
