/*! \file chkpt.cc
	\ingroup (psirb)
	\brief Ruby interface to libchkpt.
*/
#include <ruby.h>
#include "psirb.h"
#include <libpsio/psio.hpp>

extern "C" {
#include <psifiles.h>
#include <ccfiles.h>
};

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
