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
    
    ip_cwk_add(":DFT");
    
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
    
    errcod = ip_string("GRID_SPEC",&gridstring,0);
    if(errcod != IPE_OK){
	
	/* Only need one grid */
	DFT_options.grid_dim = 1;
	DFT_options.grid_info = (struct grid_info_s *) 
	    malloc(DFT_options.grid_dim*sizeof(struct grid_info_s));
	DFT_options.grid_info[0].grid_label = 
	    "Euler-McClaren / Lebedev Spheres";
	rpointstmp = 50;
	errcod = ip_data("GRID_RPOINTS","%d",&rpointstmp,0);
	DFT_options.rpoints = rpointstmp;

	angpointstmp=302;
	errcod = ip_data("GRID_ANGPOINTS","%d",&angpointstmp,0);
	DFT_options.grid_info[0].angpoints = angpointstmp;

	DFT_options.grid_info[0].leb_point = 
	    lebedev_init(DFT_options.grid_info[0].angpoints);
    }
    else
	punt("No Special Grids have been implemented in this code yet");

    /*------------------------------------
      This is for debugging purposes only
      There's no need to use separate
      print level for DFT in the future
     ------------------------------------*/
    /* Print out DFT information? */
    DFT_options.prtflag = 0;
    errcod = ip_data("DFT_PRINT","%d",&(DFT_options.prtflag),0);
    
    if(DFT_options.prtflag > 0){
	fprintf(outfile," -DFT info:\n");
	fprintf(outfile,"    Functional = %s",funcstring);
	for(i=0;i<DFT_options.grid_dim;i++){
	    fprintf(outfile,"   %s",DFT_options.grid_info[0].grid_label);
	    fprintf(outfile," Radial Cutoff = %e",
		    DFT_options.grid_info[i].rcut);
	    fprintf(outfile," Radial Points = %5d Angular Points =%5d\n",
		    DFT_options.rpoints,
		    DFT_options.grid_info[i].angpoints);
	}
	fflush(outfile);
    }
}


















