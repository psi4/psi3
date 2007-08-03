/*! \file chkpt.cc
	\ingroup (psirb)
	\brief Ruby interface to libchkpt.
*/
#include <ruby.h>
#include "psirb.h"

extern "C" {
#include <psifiles.h>
#include <ccfiles.h>
#include <libpsio/psio.hpp>
};

//! Ruby functions: Psi::Chkpt::exist? and Psi::Chkpt::exists?
/*! Ruby interface to chkpt_exist. Checks to see if the requested keyword exists in the
	checkpoint file.
	\param self Ruby object calling this function.
	\param keyword Keyword to check for.
	\return Qtrue (C) / true (Ruby) if it exists, or Qfalse (C) / false (Ruby) if it does not.
*/
VALUE ruby_psi_chkpt_exist(VALUE self, VALUE keyword)
{
	// Convert the given keyword to a C-string
	VALUE str = StringValue(keyword);
	char *p = RSTRING(str)->ptr;
	char *keyw = NULL;
	
	if (p == NULL) {
		rb_raise(rb_eArgError, "wrong argument, expected a string");
	}
	
	// Call the Psi chkpt function
	Chkpt chkpt(&Globals::g_psioDefault, PSIO_OPEN_OLD);
	keyw = chkpt.build_keyword(p);
	int result = chkpt.exist(keyw);
	free(keyw);
	
	if (result) return Qtrue;
	else        return Qfalse;
}

//! Ruby function: Psi::Chkpt::label
/*! Ruby interface to chkpt_rd_label. Reads the label from checkpoint.
	\param self Ruby object calling this function.
	\return Label as a Ruby string.
*/
VALUE ruby_psi_chkpt_label_get(VALUE self)
{
	// Read in the label from Chkpt
	char *label = NULL;
	
	Chkpt chkpt(&Globals::g_psioDefault, PSIO_OPEN_OLD);
	label = chkpt.rd_label();
	
	// Create ruby string
	VALUE result = rb_str_new2(label);
	
	return result;
}

//! Ruby function: Psi::Chkpt::escf
/*! Ruby interface to chkpt_rd_escf. Reads SCF energy from checkpoint.
	\param self Ruby object calling this function.
	\return SCF energy in a Ruby object.
*/
VALUE ruby_psi_chkpt_escf_get(VALUE self)
{
	// Read in the scf
	double escf;
	
	Chkpt chkpt(&Globals::g_psioDefault, PSIO_OPEN_OLD);
	escf = chkpt.rd_escf();
	
	VALUE result = rb_float_new(escf);
	return result;
}

//! Ruby function: Psi::Chkpt::escf=
/*! Ruby interface to chkpt_wt_escf. Writes the SCF energy to checkpoint.
	\param self Ruby object calling this function.
	\param vescf New SCF energy.
*/
VALUE ruby_psi_chkpt_escf_set(VALUE self, VALUE vescf)
{
	// Read in the scf
	double escf = NUM2DBL(vescf);
	
	Chkpt chkpt(&Globals::g_psioDefault, PSIO_OPEN_OLD);
	chkpt.wt_escf(escf);
	
	return self;
}

//! Ruby function: Psi::Chkpt::eref
/*! Ruby interface to chkpt_rd_eref. Read the reference energy from checkpoint.
	\param self Ruby object calling this function.
	\return Reference energy in a Ruby object.
*/
VALUE ruby_psi_chkpt_eref_get(VALUE self)
{
	// Read in the scf
	double escf;
	
	Chkpt chkpt(&Globals::g_psioDefault, PSIO_OPEN_OLD);
	escf = chkpt.rd_escf();
	
	VALUE result = rb_float_new(escf);
	return result;
}

//! Ruby function: Psi::Chkpt::eref=
/*! Ruby interface to chkpt_wt_ref. Writes the reference energy to checkpoint.
	\param self Ruby object calling this function.
	\param veref New reference energy as a Ruby object.
*/
VALUE ruby_psi_chkpt_eref_set(VALUE self, VALUE veref)
{
	// Read in the scf
	double eref = NUM2DBL(veref);
	
	Chkpt chkpt(&Globals::g_psioDefault, PSIO_OPEN_OLD);
	chkpt.wt_eref(eref);
	
	return self;
}

//! Ruby function: Psi::Chkpt::ecorr
/*! Ruby interface to chkpt_rd_ecorr. Read the correlation energy from checkpoint.
	\param self Ruby object calling this function.
	\return Correlation energy as a Ruby object.
*/
VALUE ruby_psi_chkpt_ecorr_get(VALUE self)
{
	double ecorr;
	
	Chkpt chkpt(&Globals::g_psioDefault, PSIO_OPEN_OLD);
	ecorr = chkpt.rd_ecorr();
	
	VALUE result = rb_float_new(ecorr);
	return result;
}

//! Ruby function: Psi::Chkpt::ecorr=
/*! Ruby interface to chkpt_wt_ecorr. Writes the new correlation energy to checkpoint.
	\param self Ruby object calling this function.
	\param vecorr New correlation energy as a Ruby object.
*/
VALUE ruby_psi_chkpt_ecorr_set(VALUE self, VALUE vecorr)
{
	double ecorr = NUM2DBL(vecorr);
	
	Chkpt chkpt(&Globals::g_psioDefault, PSIO_OPEN_OLD);
	chkpt.wt_ecorr(ecorr);
	
	return self;
}

//! Ruby function: Psi::Chkpt::enuc
/*! Ruby interface to chkpt_rd_enuc. Read the nuclear repulsion energy to checkpoint.
	\param self Ruby object calling this function.
	\return Nuclear repulsion energy as a Ruby object.
*/
VALUE ruby_psi_chkpt_enuc_get(VALUE self)
{
	double enuc;
	
	Chkpt chkpt(&Globals::g_psioDefault, PSIO_OPEN_OLD);
	enuc = chkpt.rd_enuc();
	
	VALUE result = rb_float_new(enuc);
	return result;
}

//! Ruby function: Psi::Chkpt:enuc=
/*! Ruby interface to chkpt_wt_enuc. Writes the nuclear repulsion energy to checkpoint.
	\param self Ruby object calling this function.
	\param venuc New nuclear repulsion energy as a Ruby object.
*/
VALUE ruby_psi_chkpt_enuc_set(VALUE self, VALUE venuc)
{
	double enuc = NUM2DBL(venuc);
	
	Chkpt chkpt(&Globals::g_psioDefault, PSIO_OPEN_OLD);
	chkpt.wt_enuc(enuc);
	
	return self;
}

//! Ruby function: Psi::Chkpt::efzc
/*! Ruby interface to chkpt_rd_efzc. Reads the frozen core energy from checkpoint.
	\param self Ruby object that is calling this function.
	\return Frozen core energy in a Ruby object.
*/
VALUE ruby_psi_chkpt_efzc_get(VALUE self)
{
	double efzc;
	
	Chkpt chkpt(&Globals::g_psioDefault, PSIO_OPEN_OLD);
	efzc = chkpt.rd_efzc();
	
	VALUE result = rb_float_new(efzc);
	return result;
}

//! Ruby function: Psi::Chkpt::efzc=
/*! Ruby interface to chkpt_wt_efzc. Writes the frozen core energy to checkpoint.
	\param self Ruby object that is calling this function.
	\param vefzc New frozen core energy.
*/
VALUE ruby_psi_chkpt_efzc_set(VALUE self, VALUE vefzc)
{
	double efzc = NUM2DBL(vefzc);
	
	Chkpt chkpt(&Globals::g_psioDefault, PSIO_OPEN_OLD);
	chkpt.wt_efzc(efzc);
	
	return self;
}

//! Ruby function: Psi::Chkpt::etot
/*! Ruby interface to chkpt_rd_etot. Reads the total energy from checkpoint file.
	\param self Ruby object that is calling this function.
	\return Total energy in a Ruby object.
*/
VALUE ruby_psi_chkpt_etot_get(VALUE self)
{
	double etot;
	
	Chkpt chkpt(&Globals::g_psioDefault, PSIO_OPEN_OLD);
	etot = chkpt.rd_etot();
	
	VALUE result = rb_float_new(etot);
	return result;
}

//! Ruby function: Psi::Chkpt::etot=
/*! Ruby interface to chkpt_wt_etot. Write the new total energy to checkpoint file.
	\param self Ruby object that is calling this function.
	\param vetot New total energy.
*/
VALUE ruby_psi_chkpt_etot_set(VALUE self, VALUE vetot)
{
	double etot = NUM2DBL(vetot);
	
	Chkpt chkpt(&Globals::g_psioDefault, PSIO_OPEN_OLD);
	chkpt.wt_etot(etot);
	
	return self;
}

//! Ruby function: Psi::Chkpt::disp
/*! Ruby interface to chkpt_rd_disp. Reads the current geometry displacement number.
	\param self Ruby object that is calling this function.
	\return Current displacement number.
*/
VALUE ruby_psi_chkpt_disp_get(VALUE self)
{
	int disp;
	
	Chkpt chkpt(&Globals::g_psioDefault, PSIO_OPEN_OLD);
	disp = chkpt.rd_disp();
	
	VALUE result = INT2FIX(disp);
	return result;
}

//! Ruby function: Psi::Chkpt::disp=
/*! Ruby interface to chkpt_wt_disp. Writes out the current geometry displacement number.
	\param self Ruby object that is calling this function.
	\param ndisp New displacement number.
*/
VALUE ruby_psi_chkpt_disp_set(VALUE self, VALUE ndisp)
{
	int disp = NUM2INT(ndisp);
	
	Chkpt chkpt(&Globals::g_psioDefault, PSIO_OPEN_OLD);
	chkpt.wt_disp(disp);
	
	return self;
}

//! Ruby function: Psi::Chkpt::eccsd
/*! Returns the CCSD energy contribution to the total energy.
	\param self The Ruby object that is calling this function.
	\return CCSD energy contribution as a Ruby object.
	\note Would like for this value to be saved in Chkpt and read from there.
*/
VALUE ruby_psi_chkpt_eccsd_get(VALUE self)
{
	double energy;
	
	Chkpt chkpt(&Globals::g_psioDefault, PSIO_OPEN_OLD);
	energy = chkpt.rd_eccsd();
	
	// Return the value to the user
	return rb_float_new(energy);
}

//! Ruby function: Psi::Chkpt:e_t
/*! Returns the (T) energy contribution to the total energy.
	\param self The Ruby object that is calling this function.
	\return (T) energy contribution as a Ruby object.
	\note Would like for this value to be saved in Chkpt and read from there.
*/
VALUE ruby_psi_chkpt_e_t_get(VALUE self)
{
	double energy;
	
	Chkpt chkpt(&Globals::g_psioDefault, PSIO_OPEN_OLD);
	energy = chkpt.rd_e_t();
	
	// Return the value to the user
	return rb_float_new(energy);
}

/*! Creates a module in the Psi address space for Chkpt. This address space is accessible through
	Psi::Chkpt
*/
void create_ruby_psi_chkpt_module()
{
	// Handle for Chkpt module
	VALUE rubyChkpt;
	
	// Define a sub-module of Psi named Chkpt
	rubyChkpt = rb_define_module_under(Globals::rubyPsi, "Chkpt");
	
	// Add methods to the new module
	rb_define_module_function(rubyChkpt, "exist?",  RUBYCAST(ruby_psi_chkpt_exist),     1);
	rb_define_module_function(rubyChkpt, "exists?", RUBYCAST(ruby_psi_chkpt_exist),     1);
	rb_define_module_function(rubyChkpt, "label", 	RUBYCAST(ruby_psi_chkpt_label_get), 0);
	rb_define_module_function(rubyChkpt, "escf",	RUBYCAST(ruby_psi_chkpt_escf_get),  0);
	rb_define_module_function(rubyChkpt, "escf=",	RUBYCAST(ruby_psi_chkpt_escf_set),  1);
	rb_define_module_function(rubyChkpt, "eref",	RUBYCAST(ruby_psi_chkpt_eref_get),  0);
	rb_define_module_function(rubyChkpt, "eref=", 	RUBYCAST(ruby_psi_chkpt_eref_set),  1);
	rb_define_module_function(rubyChkpt, "ecorr", 	RUBYCAST(ruby_psi_chkpt_ecorr_get), 0);
	rb_define_module_function(rubyChkpt, "ecorr=", 	RUBYCAST(ruby_psi_chkpt_ecorr_set), 1);
	rb_define_module_function(rubyChkpt, "enuc", 	RUBYCAST(ruby_psi_chkpt_enuc_get),  0);
	rb_define_module_function(rubyChkpt, "enuc=", 	RUBYCAST(ruby_psi_chkpt_enuc_set),  1);
	rb_define_module_function(rubyChkpt, "efzc", 	RUBYCAST(ruby_psi_chkpt_efzc_get),  0);
	rb_define_module_function(rubyChkpt, "efzc=", 	RUBYCAST(ruby_psi_chkpt_efzc_set),  1);
	rb_define_module_function(rubyChkpt, "etot", 	RUBYCAST(ruby_psi_chkpt_etot_get),  0);
	rb_define_module_function(rubyChkpt, "etot=", 	RUBYCAST(ruby_psi_chkpt_etot_set),  1);
	rb_define_module_function(rubyChkpt, "disp", 	RUBYCAST(ruby_psi_chkpt_etot_get),  0);
	rb_define_module_function(rubyChkpt, "disp=", 	RUBYCAST(ruby_psi_chkpt_etot_set),  1);
	rb_define_module_function(rubyChkpt, "eccsd",   RUBYCAST(ruby_psi_chkpt_eccsd_get), 0);
	rb_define_module_function(rubyChkpt, "e_t",     RUBYCAST(ruby_psi_chkpt_e_t_get),   0);
}

/////////////////////
// C++ version Chkpt
/////////////////////
Chkpt::Chkpt(psi::PSIO *psioObject, int status) : psio(psioObject)
{
	char *prefix;
	psio_tocentry *this_entry;
	
	psio->open(PSIF_CHKPT, status);

	if(psio->tocscan(PSIF_CHKPT, "Default prefix") != NULL) {
		prefix = rd_prefix();
		set_prefix(prefix);
		free(prefix);
	}
	else {
		set_prefix("");
		commit_prefix();  /* we assume that no default prefix existed in PSIF_CHKPT */
	}
}

Chkpt::~Chkpt()
{
	psio->close(PSIF_CHKPT, 1);
	psio = NULL;
}

char *Chkpt::rd_prefix()
{
	char *prefix;
	
	prefix = (char *) malloc(CHKPT_PREFIX_LEN*sizeof(char));
	
	psio->read_entry(PSIF_CHKPT, "Default prefix", prefix, CHKPT_PREFIX_LEN*sizeof(char));
	
	return prefix;
}

void Chkpt::wt_prefix(char *prefix)
{  
	psio->write_entry(PSIF_CHKPT, "Default prefix", prefix, CHKPT_PREFIX_LEN*sizeof(char));
}

void Chkpt::set_prefix(char *prefix)
{
  	strcpy(chkpt_prefix, prefix);
}

void Chkpt::commit_prefix()
{
  	wt_prefix(chkpt_prefix);
}

void Chkpt::reset_prefix()
{
	chkpt_prefix[0] = '\0';
}

char *Chkpt::get_prefix(void)
{
	char *prefix;

	prefix = (char *) malloc(CHKPT_PREFIX_LEN*sizeof(char));

	strcpy(prefix,chkpt_prefix);

	return prefix;
}

char *Chkpt::build_keyword(char *key)
{
	char *keyword;
	int keylen;

	keylen = strlen(key) + strlen(chkpt_prefix) + 2;
	if(keylen > PSIO_KEYLEN) {
		printf("LIBCHKPT: requested key exceeds allowed LIBPSIO length: :%s:%s\n", 
			chkpt_prefix, key);
		exit(PSI_RETURN_FAILURE);
	}

	keyword = (char *) malloc((keylen+1)*sizeof(char));
	sprintf(keyword, ":%s:%s", chkpt_prefix, key);
	keyword[keylen] = '\0';

	return keyword;
}

int Chkpt::exist(char *keyword)
{
	int exists=0;
	
	if (psio->tocscan(PSIF_CHKPT, keyword) != NULL)
		exists=1;
		
	return exists;
}

char *Chkpt::rd_label()
{
	char *label;
	char *keyword;
	keyword = build_keyword("Label");

	label = (char *) malloc(80 * sizeof(char));

	psio->read_entry(PSIF_CHKPT, keyword, (char *) label, 80*sizeof(char));

	free(keyword);
	return label;
}

void Chkpt::wt_label(char *label)
{
	char *keyword;
	keyword = build_keyword("Label");
	
	psio->write_entry(PSIF_CHKPT, keyword, (char*)label, 80*sizeof(char));
	
	free(keyword);
}

double Chkpt::rd_escf(void)
{
	double escf;
	char *keyword;
	keyword = build_keyword("SCF energy");

	psio->read_entry(PSIF_CHKPT, keyword, (char *) &escf, sizeof(double));

	free(keyword);
	return escf;
}

void Chkpt::wt_escf(double escf)
{
	char *keyword;
	keyword = build_keyword("SCF energy");

	psio->write_entry(PSIF_CHKPT, keyword, (char *) &escf, sizeof(double));

	free(keyword);
}

double Chkpt::rd_eref(void)
{
	double eref;
	char *keyword;
	keyword = build_keyword("Reference energy");

	psio->read_entry(PSIF_CHKPT, keyword, (char *) &eref, sizeof(double));

	free(keyword);
	return eref;
}

void Chkpt::wt_eref(double eref)
{
	char *keyword;
	keyword = build_keyword("Reference energy");

	psio->write_entry(PSIF_CHKPT, keyword, (char *) &eref, 
		sizeof(double));

	free(keyword);
}

double Chkpt::rd_ecorr(void)
{
	double ecorr;
	char *keyword;
	keyword = build_keyword("Correlation energy");

	psio->read_entry(PSIF_CHKPT, keyword, (char *) &ecorr,
		sizeof(double));

	free(keyword);
	return ecorr;
}

void Chkpt::wt_ecorr(double ecorr)
{
	char *keyword;
	keyword = build_keyword("Correlation energy");

	psio->write_entry(PSIF_CHKPT, keyword, (char *) &ecorr,
		sizeof(double));

	free(keyword);
}

double Chkpt::rd_enuc(void)
{
	double enuc;
	char *keyword;
	keyword = build_keyword("Nuclear rep. energy");

	psio->read_entry(PSIF_CHKPT, keyword, (char *) &enuc, sizeof(double));

	free(keyword);
	return enuc;
}

void Chkpt::wt_enuc(double enuc)
{
	char *keyword;
	keyword = build_keyword("Nuclear rep. energy");

	psio->write_entry(PSIF_CHKPT, keyword, (char *) &enuc, sizeof(double));

	free(keyword);
}

double Chkpt::rd_efzc(void)
{
	double efzc;
	char *keyword;
	keyword = build_keyword("Frozen core energy");

	psio->read_entry(PSIF_CHKPT, keyword, (char *) &efzc, sizeof(double));

	free(keyword);
	return efzc;
}

void Chkpt::wt_efzc(double efzc)
{
	char *keyword;
	keyword = build_keyword("Frozen core energy");

	psio->write_entry(PSIF_CHKPT, keyword, (char *) &efzc, 
		sizeof(double));

	free(keyword);
}

double Chkpt::rd_etot(void)
{
	double etot;
	char *keyword;
	keyword = build_keyword("Total energy");

	psio->read_entry(PSIF_CHKPT, keyword, (char *) &etot, sizeof(double));

	free(keyword);
	return etot;
}

void Chkpt::wt_etot(double etot)
{
	char *keyword;
	keyword = build_keyword("Total energy");

	psio->write_entry(PSIF_CHKPT, keyword, (char *) &etot, sizeof(double));

	free(keyword);
}

int Chkpt::rd_disp(void)
{
	int disp;
	char *keyword;
	keyword = build_keyword("Current displacement");

	psio->read_entry(PSIF_CHKPT, keyword, (char *) &disp,
		sizeof(int));

	free(keyword);
	return disp;
}

void Chkpt::wt_disp(int disp)
{
	char *keyword;
	keyword = build_keyword("Current displacement");

	psio->write_entry(PSIF_CHKPT, keyword, (char *) &disp,
		sizeof(int));

	free(keyword);
}

double Chkpt::rd_eccsd()
{
	double energy;
	char *keyword;
	keyword = build_keyword("CCSD Energy");
	
	// Read the energy in
	psio->read_entry(PSIF_CHKPT, keyword, (char*)&energy, sizeof(double));
		
	// Return the value to the user
	return energy;
}

void Chkpt::wt_eccsd(double eccsd)
{
	char *keyword;
	keyword = build_keyword("CCSD Energy");
	
	psio->write_entry(PSIF_CHKPT, keyword, (char*)&eccsd, sizeof(double));
	
	free(keyword);
}

double Chkpt::rd_e_t()
{
	double energy;
	char *keyword;
	keyword = build_keyword("(T) Energy");
	
	// Read the energy in
	psio->read_entry(PSIF_CHKPT, keyword, (char*)&energy, sizeof(double));
		
	// Return the value to the user
	return energy;
}

void Chkpt::wt_e_t(double e_t)
{
	char *keyword;
	keyword = build_keyword("(T) Energy");
	
	psio->write_entry(PSIF_CHKPT, keyword, (char*)&e_t, sizeof(double));
	
	free(keyword);
}
