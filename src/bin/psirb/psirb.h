#ifndef __PSIRB_H__
#define __PSIRB_H__

#ifdef MAIN
#define EXT
#else
#define EXT extern
#endif

#include <stdio.h>
#include <string>
#include <ruby.h>
#include <libpsio/psio.hpp>

#define PSI_VERSION_MAJOR 3
#define PSI_VERSION_MINOR 3

/*! Namespace for containing global variables */
namespace Globals {
/*! All output should be sent to this variable */
	EXT FILE* g_fOutput;

/*! The name of the input file, could be useful for something */
	EXT std::string g_szInputFile;
	
/*! Name of the output filename */
	EXT std::string g_szOutputFilename;
	
/*! Verbosity */
	EXT bool g_bVerbose;
		
#ifdef MAIN
/*! All classes that provide a Ruby interface need this to say which module they belong to */
	EXT VALUE g_rbPsi = Qnil;
/*! How much to indent the Ruby puts output by. */
	EXT int g_iPutsIndent = 0;
#else
	EXT VALUE g_rbPsi;
	EXT int g_iPutsIndent;
#endif
};

// Helper macros
/*! Cast a C/++ function to a Ruby acceptable form */
#define RUBYCAST(x) (VALUE (*)(...))(x)

/*! Help the user get the Psi object from Ruby. Psi object MUST have an assign function */
#define RUBYPSIDATA(RObj, OType, OPtr) \
	Data_Get_Object(RObj, OType, OPtr);

//#ifndef CHKPT_PREFIX_LEN
//#define CHKPT_PREFIX_LEN 32
//#endif

/*! A limited C++ implementation of libchkpt using C++ libpsio. */
/*
class Chkpt {
	psi::PSIO *psio;
	char chkpt_prefix[CHKPT_PREFIX_LEN];
public:
	
	Chkpt(psi::PSIO *psioObject, int status);
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
}; */

/*! A calculation will now be grouped into a Task object. Tasks provide a 
	simple way of having several prefixes in a single input file.
	There is a global Task object that is used by default. */
class Task {
	/*! The libpsio object for this task */
	psi::PSIO m_psiPSIO;
	
	/*! Equivalent to the psi_file_prefix */
	std::string m_szFilePrefix;

	/*! Default scratch space that Psi is supposed to use. */
	std::string m_szScratchPath;
	
	/*! Ruby reference to the Task class descriptor */
	static VALUE m_rbTask;
	
public:
	/*! Default constructor; sets prefix to "psi" and scratch to "/tmp/" */
	Task(std::string prefix = "psi", std::string scratch = "/tmp/");
	
	/*! Destructor */
	~Task();
	
	/*! Accessor function get for the prefix 
		\return current prefix */
	const std::string& prefix();
	/*! Accessor function put for the prefix
		\param new_prefix what to change the prefix to */
	void prefix(std::string new_prefix);
	
	/*! Accessor function get for the scratch location
		\return current location */
	const std::string& scratch();
	/*! Accessor function put for the scratch location
		\param new_scratch what to change the scratch to */
	void scratch(std::string new_scratch);
	
	//
	// Ruby framework for Task
	//
	
	/*! Creates the Ruby class framework */
	static void create_ruby_class();

	/*! Called by Ruby when it needs to delete a class. */
	static void rb_free(void *p);
	
	/*! Called by Ruby during object creation */
	static VALUE rb_alloc(VALUE klass);
	
	/*! Called by Ruby during object creation */
	static VALUE rb_init(int argc, VALUE* argv, VALUE self);
	
	/*! Called by Ruby during object copy creation */
	static VALUE rb_init_copy(VALUE copy, VALUE orig);
	
	/*! Called by Ruby if the user try to print a Task object */
	static VALUE rb_to_s(VALUE self);
	
	/*! Ruby function: Task.prefix= */
	static VALUE rb_prefix_set(VALUE self, VALUE newPrefix);
	/*! Ruby function: Task.prefix */
	static VALUE rb_prefix_get(VALUE self);

	/*! Ruby function: Task.scratch= */
	static VALUE rb_scratch_set(VALUE self, VALUE newsScratch);
	/*! Ruby function: Task.scratch */
	static VALUE rb_scratch_get(VALUE self);
	
	//
	// Checkpoint interface
	static VALUE rb_chkpt_exist(VALUE, VALUE);
	static VALUE rb_chkpt_label_get(VALUE);
	static VALUE rb_chkpt_label_set(VALUE self, VALUE label);
	static VALUE rb_chkpt_escf_get(VALUE);
	static VALUE rb_chkpt_escf_set(VALUE, VALUE);
	static VALUE rb_chkpt_eref_get(VALUE);
	static VALUE rb_chkpt_eref_set(VALUE, VALUE);
	static VALUE rb_chkpt_ecorr_get(VALUE);
	static VALUE rb_chkpt_ecorr_set(VALUE, VALUE);
	static VALUE rb_chkpt_enuc_get(VALUE);
	static VALUE rb_chkpt_enuc_set(VALUE, VALUE);
	static VALUE rb_chkpt_efzc_get(VALUE);
	static VALUE rb_chkpt_efzc_set(VALUE, VALUE);
	static VALUE rb_chkpt_etot_get(VALUE);
	static VALUE rb_chkpt_etot_set(VALUE, VALUE);
	static VALUE rb_chkpt_disp_get(VALUE);
	static VALUE rb_chkpt_disp_set(VALUE, VALUE);
	static VALUE rb_chkpt_eccsd_get(VALUE);
	static VALUE rb_chkpt_e_t_get(VALUE);
	static VALUE rb_chkpt_emp2_get(VALUE);	
};

/*
class ZEntry {
private:
	z_entry *m_zEntry;

	//! Ruby reference to the Z-Matrix class descriptor
	static VALUE m_rbZEntry;
	
public:
	//! Default constructor
	ZEntry();
	~ZEntry();
	
	//
	// Ruby framework for ZEntry
	//
	
	//! Creates the Ruby class framework
	static void create_ruby_class();
	
	//! Called by Ruby when it needs to delete a class.
	
};
*/

#endif // __PSIRB_H__
