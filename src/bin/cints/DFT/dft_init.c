#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include"functional.h"
#include"calc_den.h"
#include"lebedev_init.h"
#include"grid_init.h"

void dft_init(void){
    
    int errcod;
    int i,j;
    int si,sj; 
    int si_n_ao;
    int sj_n_ao;
    int rpointstmp;
    int angpointstmp;
    int n_shells;
    int elem_in_func;
    char *exchstring;
    char *corrstring;
    char *gridstring;
    char *funcstring;
    double tmpa, tmpb;
    struct shell_pair *sp;
    
    /* DFT input parsing */
    elem_in_func = 0;
    errcod = ip_count("FUNCTIONAL",&elem_in_func,0);
    if(elem_in_func == 0){
	errcod=ip_string("FUNCTIONAL",&funcstring,0);
	if(errcod == IPE_OK){
	    if(!strcmp(funcstring,"XALPHA")){
		UserOptions.hf_exch = 0.0;
		DFT_options.exchange_function = slater;
		DFT_options.exchange_V_function = d_slater;
		DFT_options.den_calc = calc_density;
		DFT_options.correlation_function = no_funct;
		DFT_options.correlation_V_function = no_funct;
	    }
	    else
		punt("Unrecognized fucntional specified with keyword FUNCTIONAL");
	}
	else
	    punt("Must define a functional with keyword FUNCTIONAL");
    }
    
    else if(elem_in_func == 2){
	errcod = ip_string("FUNCTIONAL",&exchstring,1,0);
	
        if(!strcmp(exchstring,"SLATER")){
	    UserOptions.hf_exch = 0.0;
	    DFT_options.exchange_function = slater;
	    DFT_options.exchange_V_function = d_slater;
	    DFT_options.den_calc = calc_density;
	}
	else if(!strcmp(exchstring,"DENSITY")){
	    DFT_options.exchange_function = density;
	    DFT_options.exchange_V_function = density;
	    DFT_options.den_calc = calc_density;
	}
	else
	    punt("Unrecognized or nonimplemented exchange functional specified");
	
	errcod = ip_string("FUNCTIONAL",&corrstring,1,1);
	
	if(!strcmp(corrstring,"NONE")){
	    DFT_options.correlation_function = no_funct;
	    DFT_options.correlation_V_function = no_funct;
	    corrstring[0] = '\0';
	}
        else if(!strcmp(corrstring,"VWN5")){
            DFT_options.correlation_function = VWN5_r;
            DFT_options.correlation_V_function = d_VWN5_r;
        }           
	else if(!strcmp(corrstring,"VWN4")){
            DFT_options.correlation_function = VWN4_r;
            DFT_options.correlation_V_function = d_VWN4_r;
        }           
	else
	    punt("Unrecognized or nonimplemented correlation fuctional specified");

	i = strlen(exchstring) + strlen(corrstring);
	if (i == 0)
	  funcstring = strdup("NONE");
	else {
	  funcstring = (char *) malloc(sizeof(char)* (i + 1));
	  sprintf(funcstring,"%s%s",exchstring,corrstring);
	}
	free(exchstring); free(corrstring);
    }
    else
	punt("Something wrong in the specification of the FUNCTIONAL keyword");    
}

void cleanup_dft_options(DFT_options_t DFT_options){
    
    free(DFT_options.basis);
    free(DFT_options.Bragg);
    cleanup_grid_type(DFT_options.grid);
}
















