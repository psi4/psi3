/*! \file scf_parsing.c
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include<stdio.h>
#include<math.h>
#include<libipv1/ip_lib.h>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include"dft_init.h"

void scf_parsing(void)
{
    int errcod;
    char *refstring;

    UserOptions.make_oei = 0;
    UserOptions.make_eri = 0;
    UserOptions.make_fock = 1;
    UserOptions.print_lvl = 0;
    UserOptions.symm_ints = 0;
    UserOptions.make_dft = 0;
    UserOptions.hf_exch = 1.0;
    errcod = ip_string("REFERENCE",&refstring,0);
    if (errcod != IPE_OK)
	punt("REFERENCE keyword is missing");
    else if (!strcmp(refstring,"RHF") || !strcmp(refstring,""))
	UserOptions.reftype = rhf;
    else if (!strcmp(refstring,"ROHF"))
	UserOptions.reftype = rohf;
    else if (!strcmp(refstring,"UHF"))
	UserOptions.reftype = uhf;
    else if (!strcmp(refstring,"RKS")){
	UserOptions.reftype = rhf;
	UserOptions.make_dft = 1;
	dft_init();
    }
    else if (!strcmp(refstring,"UKS")){
	UserOptions.reftype = uhf;
	UserOptions.make_dft = 1;
	dft_init();
    }
    else
	punt("The specified REFERENCE not implemented");
}
