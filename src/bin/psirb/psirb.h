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
	
/*! Equivalent to the psi_file_prefix */
	EXT std::string g_szFilePrefix;
	
/*! Default scratch space that Psi is supposed to use. */
	EXT std::string g_szScratchPath;
	
/*! Default psio object. This is for simple input files. */
	EXT psi::PSIO g_psioDefault;

/*! Mapping a PSIO object to a string name. Ruby passes a 
	to the Ruby C-functions. The C-functions use this
	mapping to get the right PSIO object. */
	EXT std::map<std::string, psi::PSIO*> g_mapPSIO;
	
#ifdef MAIN
/*! All classes that provide a Ruby interface need this to say which module they belong to */
	EXT VALUE rubyPsi = Qnil;
/*! How much to indent the Ruby puts output by. */
	EXT int g_iPutsIndent = 0;
#else
	EXT VALUE rubyPsi;
	EXT int g_iPutsIndent;
#endif
};

// Helper macros
/*! Cast a C/++ function to a Ruby acceptable form */
#define RUBYCAST(x) (VALUE (*)(...))(x)

/*! Help the user get the Psi object from Ruby. Psi object MUST have an assign function */
#define RUBYPSIDATA(RObj, OType, OPtr) \
	Data_Get_Object(RObj, OType, OPtr);

#ifndef CHKPT_PREFIX_LEN
#define CHKPT_PREFIX_LEN 32
#endif

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
};

#endif // __PSIRB_H__
