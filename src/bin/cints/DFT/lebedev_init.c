/* ----------------------------------------
   
   lebedev_init.c
   
   by Shawn Brown
   
   This will initialize the values for the
   lebedev grid

   ---------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <memory.h>
#include <stdlib.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <libint.h>
#include <pthread.h>

#include"defines.h"
#define EXTERN
#include"global.h"

struct leb_point_s *lebedev_init(int degree){
    
    int i,j,ij;
    int offset;
    double *As,*Bs,*Cs,*Ds;
    struct leb_point_s *leb_tmp;
    
    leb_tmp = (struct leb_point_s *)
	malloc(degree*sizeof(struct leb_point_s));
    
    switch (degree){
	
    case 302:
	
/* a1 points */
	As = init_array(2);
	Bs = init_array(6);
	Cs = init_array(2);
	Ds = init_array(2);
	
	As[0] = 8.54591172878e-4;
	As[1] = 0.0;
	As[2] = 3.59911928502e-3;
	
	Bs[0] = 3.65004580768e-3;
	Bs[1] = 3.60482260142e-3;
	Bs[2] = 3.57672966173e-3;
	Bs[3] = 3.44978842429e-3;
	Bs[4] = 3.10895312238e-3;
	Bs[5] = 2.35210141366e-3;
	
	Cs[0] = 3.60082093222e-3;
	Cs[1] = 2.98234496317e-3;
	
	Ds[0] = 3.57154055427e-3;
	Ds[1] = 3.39231220501e-3;
	
	for(i=0;i<6;i++)
	    leb_tmp[i].ang_quad_weight = As[0];
	for(i=6;i<14;i++)
	    leb_tmp[i].ang_quad_weight = As[2];
	offset = 14;
	for(i=0;i<6;i++){
	    for(j=0;j<24;j++){
		ij = ((i*24)+j)+offset;
		leb_tmp[ij].ang_quad_weight = Bs[i];
	    }
	}
	offset = 158;
	for(i=0;i<2;i++){
	    for(j=0;j<24;j++){
		ij = ((i*24)+j)+offset;
		leb_tmp[ij].ang_quad_weight = Cs[i];
	    }
	}
	offset = 206;
	for(i=0;i<2;i++){
	    for(j=0;j<48;j++){
		ij = ((i*48)+j)+offset;
		leb_tmp[ij].ang_quad_weight = Ds[i];
	    }
	}
	
	leb_tmp[0].p_cart.x = 1.0;
	leb_tmp[0].p_cart.y = 0.0;
	leb_tmp[0].p_cart.z = 0.0;
	
	leb_tmp[1].p_cart.x = 0.0;
	leb_tmp[1].p_cart.y = 1.0;
	leb_tmp[1].p_cart.z = 0.0;
	
	leb_tmp[2].p_cart.x = 0.0;
	leb_tmp[2].p_cart.y = 0.0;
	leb_tmp[2].p_cart.z = 1.0;
	
	leb_tmp[3].p_cart.x = -1.0;
	leb_tmp[3].p_cart.y = 0.0;
	leb_tmp[3].p_cart.z = 0.0;
	
	leb_tmp[4].p_cart.x = 0.0;
	leb_tmp[4].p_cart.y = -1.0;
	leb_tmp[4].p_cart.z = 0.0;
	
	leb_tmp[5].p_cart.x = 0.0;
	leb_tmp[5].p_cart.y = 0.0;
	leb_tmp[5].p_cart.z = -1.0;                                            
      
	/* a3 points */
	
	leb_tmp[6].p_cart.x = -0.577350269189625842;
	leb_tmp[6].p_cart.y = 0.577350269189625842;
	leb_tmp[6].p_cart.z = 0.577350269189625842;
	
	leb_tmp[7].p_cart.x = 0.577350269189625842;
	leb_tmp[7].p_cart.y = -0.577350269189625842;
	leb_tmp[7].p_cart.z = 0.577350269189625842;
	
	leb_tmp[8].p_cart.x = 0.577350269189625842;
	leb_tmp[8].p_cart.y = 0.577350269189625842;
	leb_tmp[8].p_cart.z = -0.577350269189625842; 
	
	leb_tmp[9].p_cart.x = -0.577350269189625842;
	leb_tmp[9].p_cart.y = -0.577350269189625842;
	leb_tmp[9].p_cart.z = 0.577350269189625842; 
	
	leb_tmp[10].p_cart.x = -0.577350269189625842;
	leb_tmp[10].p_cart.y = 0.577350269189625842;
	leb_tmp[10].p_cart.z = -0.577350269189625842; 
	
	leb_tmp[11].p_cart.x = 0.577350269189625842;
	leb_tmp[11].p_cart.y = -0.577350269189625842;
	leb_tmp[11].p_cart.z = -0.577350269189625842; 
	
	leb_tmp[12].p_cart.x = -0.577350269189625842;
	leb_tmp[12].p_cart.y = -0.577350269189625842;
	leb_tmp[12].p_cart.z = -0.577350269189625842; 
	
	leb_tmp[13].p_cart.x = 0.577350269189625842;
	leb_tmp[13].p_cart.y = 0.577350269189625842;
	leb_tmp[13].p_cart.z = 0.577350269189625842; 
	
	/* b11 points */
	leb_tmp[14].p_cart.x = 0.701176641609;
	leb_tmp[14].p_cart.y = 0.701176641609;
	leb_tmp[14].p_cart.z = 0.129238672710;
	
	leb_tmp[15].p_cart.x = -0.701176641609;
	leb_tmp[15].p_cart.y = 0.701176641609;
	leb_tmp[15].p_cart.z = 0.129238672710;
	
	leb_tmp[16].p_cart.x = 0.701176641609;
	leb_tmp[16].p_cart.y = -0.701176641609;
	leb_tmp[16].p_cart.z = 0.129238672710;
 
	leb_tmp[17].p_cart.x = 0.701176641609;
	leb_tmp[17].p_cart.y = 0.701176641609;
	leb_tmp[17].p_cart.z = -0.129238672710;
	
	leb_tmp[18].p_cart.x = -0.701176641609;
	leb_tmp[18].p_cart.y = -0.701176641609;
	leb_tmp[18].p_cart.z = 0.129238672710;
	
	leb_tmp[19].p_cart.x = -0.701176641609;
	leb_tmp[19].p_cart.y = 0.701176641609;
	leb_tmp[19].p_cart.z = -0.129238672710;
	
	leb_tmp[20].p_cart.x = 0.701176641609;
	leb_tmp[20].p_cart.y = -0.701176641609;
	leb_tmp[20].p_cart.z = -0.129238672710;
	
	leb_tmp[21].p_cart.x = -0.701176641609;
	leb_tmp[21].p_cart.y = -0.701176641609;
	leb_tmp[21].p_cart.z = -0.129238672710;
	
	leb_tmp[22].p_cart.x = 0.701176641609;
	leb_tmp[22].p_cart.y = 0.129238672710;
	leb_tmp[22].p_cart.z = 0.701176641609;
	
	leb_tmp[23].p_cart.x = -0.701176641609;
	leb_tmp[23].p_cart.y = 0.129238672710;
	leb_tmp[23].p_cart.z = 0.701176641609;
	
	leb_tmp[24].p_cart.x = 0.701176641609;
	leb_tmp[24].p_cart.y = -0.129238672710;
	leb_tmp[24].p_cart.z = 0.701176641609;
	
	leb_tmp[25].p_cart.x = 0.701176641609;
	leb_tmp[25].p_cart.y = 0.129238672710;
	leb_tmp[25].p_cart.z = -0.701176641609;
	
	leb_tmp[26].p_cart.x = -0.701176641609;
	leb_tmp[26].p_cart.y = -0.129238672710;
	leb_tmp[26].p_cart.z = 0.701176641609;
	
	leb_tmp[27].p_cart.x = -0.701176641609;
	leb_tmp[27].p_cart.y = 0.129238672710;
	leb_tmp[27].p_cart.z = -0.701176641609;
	
	leb_tmp[28].p_cart.x = 0.701176641609;
	leb_tmp[28].p_cart.y = -0.129238672710;
	leb_tmp[28].p_cart.z = -0.701176641609;
	
	leb_tmp[29].p_cart.x = -0.701176641609;
	leb_tmp[29].p_cart.y = -0.129238672710;
	leb_tmp[29].p_cart.z = -0.701176641609;
	
	leb_tmp[30].p_cart.x = 0.129238672710;
	leb_tmp[30].p_cart.y = 0.701176641609;
	leb_tmp[30].p_cart.z = 0.701176641609;

	leb_tmp[31].p_cart.x = -0.129238672710;
	leb_tmp[31].p_cart.y = 0.701176641609;
	leb_tmp[31].p_cart.z = 0.701176641609;

	leb_tmp[32].p_cart.x = 0.129238672710;
	leb_tmp[32].p_cart.y = -0.701176641609;
	leb_tmp[32].p_cart.z = 0.701176641609;

	leb_tmp[33].p_cart.x = 0.129238672710;
	leb_tmp[33].p_cart.y = 0.701176641609;
	leb_tmp[33].p_cart.z = -0.701176641609;

	leb_tmp[34].p_cart.x = -0.129238672710;
	leb_tmp[34].p_cart.y = -0.701176641609;
	leb_tmp[34].p_cart.z = 0.701176641609;

	leb_tmp[35].p_cart.x = -0.129238672710;
	leb_tmp[35].p_cart.y = 0.701176641609;
	leb_tmp[35].p_cart.z = -0.701176641609;

	leb_tmp[36].p_cart.x = 0.129238672710;
	leb_tmp[36].p_cart.y = -0.701176641609;
	leb_tmp[36].p_cart.z = -0.701176641609;

	leb_tmp[37].p_cart.x = -0.129238672710;
	leb_tmp[37].p_cart.y = -0.701176641609;
	leb_tmp[37].p_cart.z = -0.701176641609;

	leb_tmp[38].p_cart.x = 0.656632941022;
	leb_tmp[38].p_cart.y = 0.656632941022;
	leb_tmp[38].p_cart.z = 0.371034178385;

	leb_tmp[39].p_cart.x = -0.656632941022;
	leb_tmp[39].p_cart.y = 0.656632941022;
	leb_tmp[39].p_cart.z = 0.371034178385;

	leb_tmp[40].p_cart.x = 0.656632941022;
	leb_tmp[40].p_cart.y = -0.656632941022;
	leb_tmp[40].p_cart.z = 0.371034178385;

	leb_tmp[41].p_cart.x = 0.656632941022;
	leb_tmp[41].p_cart.y = 0.656632941022;
	leb_tmp[41].p_cart.z = -0.371034178385;

	leb_tmp[42].p_cart.x = -0.656632941022;
	leb_tmp[42].p_cart.y = -0.656632941022;
	leb_tmp[42].p_cart.z = 0.371034178385;

	leb_tmp[43].p_cart.x = -0.656632941022;
	leb_tmp[43].p_cart.y = 0.656632941022;
	leb_tmp[43].p_cart.z = -0.371034178385;

	leb_tmp[44].p_cart.x = 0.656632941022;
	leb_tmp[44].p_cart.y = -0.656632941022;
	leb_tmp[44].p_cart.z = -0.371034178385;

	leb_tmp[45].p_cart.x = -0.656632941022;
	leb_tmp[45].p_cart.y = -0.656632941022;
	leb_tmp[45].p_cart.z = -0.371034178385;

	leb_tmp[46].p_cart.x = 0.656632941022;
	leb_tmp[46].p_cart.y = 0.371034178385;
	leb_tmp[46].p_cart.z = 0.656632941022;

	leb_tmp[47].p_cart.x = -0.656632941022;
	leb_tmp[47].p_cart.y = 0.371034178385;
	leb_tmp[47].p_cart.z = 0.656632941022;

	leb_tmp[48].p_cart.x = 0.656632941022;
	leb_tmp[48].p_cart.y = -0.371034178385;
	leb_tmp[48].p_cart.z = 0.656632941022;

	leb_tmp[49].p_cart.x = 0.656632941022;
	leb_tmp[49].p_cart.y = 0.371034178385;
	leb_tmp[49].p_cart.z = -0.656632941022;

	leb_tmp[50].p_cart.x = -0.656632941022;
	leb_tmp[50].p_cart.y = -0.371034178385;
	leb_tmp[50].p_cart.z = 0.656632941022;

	leb_tmp[51].p_cart.x = -0.656632941022;
	leb_tmp[51].p_cart.y = 0.371034178385;
	leb_tmp[51].p_cart.z = -0.656632941022;

	leb_tmp[52].p_cart.x = 0.656632941022;
	leb_tmp[52].p_cart.y = -0.371034178385;
	leb_tmp[52].p_cart.z = -0.656632941022;

	leb_tmp[53].p_cart.x = -0.656632941022;
	leb_tmp[53].p_cart.y = -0.371034178385;
	leb_tmp[53].p_cart.z = -0.656632941022;

	leb_tmp[54].p_cart.x = 0.371034178385;
	leb_tmp[54].p_cart.y = 0.656632941022;
	leb_tmp[54].p_cart.z = 0.656632941022;

	leb_tmp[55].p_cart.x = -0.371034178385;
	leb_tmp[55].p_cart.y = 0.656632941022;
	leb_tmp[55].p_cart.z = 0.656632941022;

	leb_tmp[56].p_cart.x = 0.371034178385;
	leb_tmp[56].p_cart.y = -0.656632941022;
	leb_tmp[56].p_cart.z = 0.656632941022;

	leb_tmp[57].p_cart.x = 0.371034178385;
	leb_tmp[57].p_cart.y = 0.656632941022;
	leb_tmp[57].p_cart.z = -0.656632941022;

	leb_tmp[58].p_cart.x = -0.371034178385;
	leb_tmp[58].p_cart.y = -0.656632941022;
	leb_tmp[58].p_cart.z = 0.656632941022;

	leb_tmp[59].p_cart.x = -0.371034178385;
	leb_tmp[59].p_cart.y = 0.656632941022;
	leb_tmp[59].p_cart.z = -0.656632941022;

	leb_tmp[60].p_cart.x = 0.371034178385;
	leb_tmp[60].p_cart.y = -0.656632941022;
	leb_tmp[60].p_cart.z = -0.656632941022;

	leb_tmp[61].p_cart.x = -0.371034178385;
	leb_tmp[61].p_cart.y = -0.656632941022;
	leb_tmp[61].p_cart.z = -0.656632941022;

	leb_tmp[62].p_cart.x = 0.472905413258;
	leb_tmp[62].p_cart.y = 0.472905413258;
	leb_tmp[62].p_cart.z = 0.743452042987;

	leb_tmp[63].p_cart.x = -0.472905413258;
	leb_tmp[63].p_cart.y = 0.472905413258;
	leb_tmp[63].p_cart.z = 0.743452042987;

	leb_tmp[64].p_cart.x = 0.472905413258;
	leb_tmp[64].p_cart.y = -0.472905413258;
	leb_tmp[64].p_cart.z = 0.743452042987;

	leb_tmp[65].p_cart.x = 0.472905413258;
	leb_tmp[65].p_cart.y = 0.472905413258;
	leb_tmp[65].p_cart.z = -0.743452042987;

	leb_tmp[66].p_cart.x = -0.472905413258;
	leb_tmp[66].p_cart.y = -0.472905413258;
	leb_tmp[66].p_cart.z = 0.743452042987;

	leb_tmp[67].p_cart.x = -0.472905413258;
	leb_tmp[67].p_cart.y = 0.472905413258;
	leb_tmp[67].p_cart.z = -0.743452042987;

	leb_tmp[68].p_cart.x = 0.472905413258;
	leb_tmp[68].p_cart.y = -0.472905413258;
	leb_tmp[68].p_cart.z = -0.743452042987;

	leb_tmp[69].p_cart.x = -0.472905413258;
	leb_tmp[69].p_cart.y = -0.472905413258;
	leb_tmp[69].p_cart.z = -0.743452042987;

	leb_tmp[70].p_cart.x = 0.472905413258;
	leb_tmp[70].p_cart.y = 0.743452042987;
	leb_tmp[70].p_cart.z = 0.472905413258;

	leb_tmp[71].p_cart.x = -0.472905413258;
	leb_tmp[71].p_cart.y = 0.743452042987;
	leb_tmp[71].p_cart.z = 0.472905413258;

	leb_tmp[72].p_cart.x = 0.472905413258;
	leb_tmp[72].p_cart.y = -0.743452042987;
	leb_tmp[72].p_cart.z = 0.472905413258;

	leb_tmp[73].p_cart.x = 0.472905413258;
	leb_tmp[73].p_cart.y = 0.743452042987;
	leb_tmp[73].p_cart.z = -0.472905413258;

	leb_tmp[74].p_cart.x = -0.472905413258;
	leb_tmp[74].p_cart.y = -0.743452042987;
	leb_tmp[74].p_cart.z = 0.472905413258;

	leb_tmp[75].p_cart.x = -0.472905413258;
	leb_tmp[75].p_cart.y = 0.743452042987;
	leb_tmp[75].p_cart.z = -0.472905413258;
	
	leb_tmp[76].p_cart.x = 0.472905413258;
	leb_tmp[76].p_cart.y = -0.743452042987;
	leb_tmp[76].p_cart.z = -0.472905413258;
	
	leb_tmp[77].p_cart.x = -0.472905413258;
	leb_tmp[77].p_cart.y = -0.743452042987;
	leb_tmp[77].p_cart.z = -0.472905413258;
	
	leb_tmp[78].p_cart.x = 0.743452042987;
	leb_tmp[78].p_cart.y = 0.472905413258;
	leb_tmp[78].p_cart.z = 0.472905413258;
	
	leb_tmp[79].p_cart.x = -0.743452042987;
	leb_tmp[79].p_cart.y = 0.472905413258;
	leb_tmp[79].p_cart.z = 0.472905413258;
	
	leb_tmp[80].p_cart.x = 0.743452042987;
	leb_tmp[80].p_cart.y = -0.472905413258;
	leb_tmp[80].p_cart.z = 0.472905413258;
	
	leb_tmp[81].p_cart.x = 0.743452042987;
	leb_tmp[81].p_cart.y = 0.472905413258;
	leb_tmp[81].p_cart.z = -0.472905413258;
	
	leb_tmp[82].p_cart.x = -0.743452042987;
	leb_tmp[82].p_cart.y = -0.472905413258;
	leb_tmp[82].p_cart.z = 0.472905413258;
	
	leb_tmp[83].p_cart.x = -0.743452042987;
	leb_tmp[83].p_cart.y = 0.472905413258;
	leb_tmp[83].p_cart.z = -0.472905413258;
	
	leb_tmp[84].p_cart.x = 0.743452042987;
	leb_tmp[84].p_cart.y = -0.472905413258;
	leb_tmp[84].p_cart.z = -0.472905413258;
	
	leb_tmp[85].p_cart.x = -0.743452042987;
	leb_tmp[85].p_cart.y = -0.472905413258;
	leb_tmp[85].p_cart.z = -0.472905413258;
	
	leb_tmp[86].p_cart.x = 0.351564034558;
	leb_tmp[86].p_cart.y = 0.351564034558;
	leb_tmp[86].p_cart.z = 0.867643624544;

	leb_tmp[87].p_cart.x = -0.351564034558;
	leb_tmp[87].p_cart.y = 0.351564034558;
	leb_tmp[87].p_cart.z = 0.867643624544;

	leb_tmp[88].p_cart.x = 0.351564034558;
	leb_tmp[88].p_cart.y = -0.351564034558;
	leb_tmp[88].p_cart.z = 0.867643624544;

	leb_tmp[89].p_cart.x = 0.351564034558;
	leb_tmp[89].p_cart.y = 0.351564034558;
	leb_tmp[89].p_cart.z = -0.867643624544;

	leb_tmp[90].p_cart.x = -0.351564034558;
	leb_tmp[90].p_cart.y = -0.351564034558;
	leb_tmp[90].p_cart.z = 0.867643624544;

	leb_tmp[91].p_cart.x = -0.351564034558;
	leb_tmp[91].p_cart.y = 0.351564034558;
	leb_tmp[91].p_cart.z = -0.867643624544;

	leb_tmp[92].p_cart.x = 0.351564034558;
	leb_tmp[92].p_cart.y = -0.351564034558;
	leb_tmp[92].p_cart.z = -0.867643624544;

	leb_tmp[93].p_cart.x = -0.351564034558;
	leb_tmp[93].p_cart.y = -0.351564034558;
	leb_tmp[93].p_cart.z = -0.867643624544;

	leb_tmp[94].p_cart.x = 0.351564034558;
	leb_tmp[94].p_cart.y = 0.867643624544;
	leb_tmp[94].p_cart.z = 0.351564034558;

	leb_tmp[95].p_cart.x = -0.351564034558;
	leb_tmp[95].p_cart.y = 0.867643624544;
	leb_tmp[95].p_cart.z = 0.351564034558;

	leb_tmp[96].p_cart.x = 0.351564034558;
	leb_tmp[96].p_cart.y = -0.867643624544;
	leb_tmp[96].p_cart.z = 0.351564034558;

	leb_tmp[97].p_cart.x = 0.351564034558;
	leb_tmp[97].p_cart.y = 0.867643624544;
	leb_tmp[97].p_cart.z = -0.351564034558;

	leb_tmp[98].p_cart.x = -0.351564034558;
	leb_tmp[98].p_cart.y = -0.867643624544;
	leb_tmp[98].p_cart.z = 0.351564034558;

	leb_tmp[99].p_cart.x = -0.351564034558;
	leb_tmp[99].p_cart.y = 0.867643624544;
	leb_tmp[99].p_cart.z = -0.351564034558;

	leb_tmp[100].p_cart.x = 0.351564034558;
	leb_tmp[100].p_cart.y = -0.867643624544;
	leb_tmp[100].p_cart.z = -0.351564034558;

	leb_tmp[101].p_cart.x = -0.351564034558;
	leb_tmp[101].p_cart.y = -0.867643624544;
	leb_tmp[101].p_cart.z = -0.351564034558;

	leb_tmp[102].p_cart.x = 0.867643624544;
	leb_tmp[102].p_cart.y = 0.351564034558;
	leb_tmp[102].p_cart.z = 0.351564034558;

	leb_tmp[103].p_cart.x = -0.867643624544;
	leb_tmp[103].p_cart.y = 0.351564034558;
	leb_tmp[103].p_cart.z = 0.351564034558;

	leb_tmp[104].p_cart.x = 0.867643624544;
	leb_tmp[104].p_cart.y = -0.351564034558;
	leb_tmp[104].p_cart.z = 0.351564034558;

	leb_tmp[105].p_cart.x = 0.867643624544;
	leb_tmp[105].p_cart.y = 0.351564034558;
	leb_tmp[105].p_cart.z = -0.351564034558;

	leb_tmp[106].p_cart.x = -0.867643624544;
	leb_tmp[106].p_cart.y = -0.351564034558;
	leb_tmp[106].p_cart.z = 0.351564034558;

	leb_tmp[107].p_cart.x = -0.867643624544;
	leb_tmp[107].p_cart.y = 0.351564034558;
	leb_tmp[107].p_cart.z = -0.351564034558;

	leb_tmp[108].p_cart.x = 0.867643624544;
	leb_tmp[108].p_cart.y = -0.351564034558;
	leb_tmp[108].p_cart.z = -0.351564034558;

	leb_tmp[109].p_cart.x = -0.867643624544;
	leb_tmp[109].p_cart.y = -0.351564034558;
	leb_tmp[109].p_cart.z = -0.351564034558;

	leb_tmp[110].p_cart.x = 0.221964523631;
	leb_tmp[110].p_cart.y = 0.221964523631;
	leb_tmp[110].p_cart.z = 0.949454317226;

	leb_tmp[111].p_cart.x = -0.221964523631;
	leb_tmp[111].p_cart.y = 0.221964523631;
	leb_tmp[111].p_cart.z = 0.949454317226;

	leb_tmp[112].p_cart.x = 0.221964523631;
	leb_tmp[112].p_cart.y = -0.221964523631;
	leb_tmp[112].p_cart.z = 0.949454317226;

	leb_tmp[113].p_cart.x = 0.221964523631;
	leb_tmp[113].p_cart.y = 0.221964523631;
	leb_tmp[113].p_cart.z = -0.949454317226;

	leb_tmp[114].p_cart.x = -0.221964523631;
	leb_tmp[114].p_cart.y = -0.221964523631;
	leb_tmp[114].p_cart.z = 0.949454317226;

	leb_tmp[115].p_cart.x = -0.221964523631;
	leb_tmp[115].p_cart.y = 0.221964523631;
	leb_tmp[115].p_cart.z = -0.949454317226;

	leb_tmp[116].p_cart.x = 0.221964523631;
	leb_tmp[116].p_cart.y = -0.221964523631;
	leb_tmp[116].p_cart.z = -0.949454317226;

	leb_tmp[117].p_cart.x = -0.221964523631;
	leb_tmp[117].p_cart.y = -0.221964523631;
	leb_tmp[117].p_cart.z = -0.949454317226;

	leb_tmp[118].p_cart.x = 0.221964523631;
	leb_tmp[118].p_cart.y = 0.949454317226;
	leb_tmp[118].p_cart.z = 0.221964523631;

	leb_tmp[119].p_cart.x = -0.221964523631;
	leb_tmp[119].p_cart.y = 0.949454317226;
	leb_tmp[119].p_cart.z = 0.221964523631;

	leb_tmp[120].p_cart.x = 0.221964523631;
	leb_tmp[120].p_cart.y = -0.949454317226;
	leb_tmp[120].p_cart.z = 0.221964523631;

	leb_tmp[121].p_cart.x = 0.221964523631;
	leb_tmp[121].p_cart.y = 0.949454317226;
	leb_tmp[121].p_cart.z = -0.221964523631;

	leb_tmp[122].p_cart.x = -0.221964523631;
	leb_tmp[122].p_cart.y = -0.949454317226;
	leb_tmp[122].p_cart.z = 0.221964523631;

	leb_tmp[123].p_cart.x = -0.221964523631;
	leb_tmp[123].p_cart.y = 0.949454317226;
	leb_tmp[123].p_cart.z = -0.221964523631;

	leb_tmp[124].p_cart.x = 0.221964523631;
	leb_tmp[124].p_cart.y = -0.949454317226;
	leb_tmp[124].p_cart.z = -0.221964523631;

	leb_tmp[125].p_cart.x = -0.221964523631;
	leb_tmp[125].p_cart.y = -0.949454317226;
	leb_tmp[125].p_cart.z = -0.221964523631;

	leb_tmp[126].p_cart.x = 0.949454317226;
	leb_tmp[126].p_cart.y = 0.221964523631;
	leb_tmp[126].p_cart.z = 0.221964523631;

	leb_tmp[127].p_cart.x = -0.949454317226;
	leb_tmp[127].p_cart.y = 0.221964523631;
	leb_tmp[127].p_cart.z = 0.221964523631;

	leb_tmp[128].p_cart.x = 0.949454317226;
	leb_tmp[128].p_cart.y = -0.221964523631;
	leb_tmp[128].p_cart.z = 0.221964523631;

	leb_tmp[129].p_cart.x = 0.949454317226;
	leb_tmp[129].p_cart.y = 0.221964523631;
	leb_tmp[129].p_cart.z = -0.221964523631;

	leb_tmp[130].p_cart.x = -0.949454317226;
	leb_tmp[130].p_cart.y = -0.221964523631;
	leb_tmp[130].p_cart.z = 0.221964523631;

	leb_tmp[131].p_cart.x = -0.949454317226;
	leb_tmp[131].p_cart.y = 0.221964523631;
	leb_tmp[131].p_cart.z = -0.221964523631;

	leb_tmp[132].p_cart.x = 0.949454317226;
	leb_tmp[132].p_cart.y = -0.221964523631;
	leb_tmp[132].p_cart.z = -0.221964523631;

	leb_tmp[133].p_cart.x = -0.949454317226;
	leb_tmp[133].p_cart.y = -0.221964523631;
	leb_tmp[133].p_cart.z = -0.221964523631;

	leb_tmp[134].p_cart.x = 0.096183085230;
	leb_tmp[134].p_cart.y = 0.096183085230;
	leb_tmp[134].p_cart.z = 0.990705621379;

	leb_tmp[135].p_cart.x = -0.096183085230;
	leb_tmp[135].p_cart.y = 0.096183085230;
	leb_tmp[135].p_cart.z = 0.990705621379;

	leb_tmp[136].p_cart.x = 0.096183085230;
	leb_tmp[136].p_cart.y = -0.096183085230;
	leb_tmp[136].p_cart.z = 0.990705621379;

	leb_tmp[137].p_cart.x = 0.096183085230;
	leb_tmp[137].p_cart.y = 0.096183085230;
	leb_tmp[137].p_cart.z = -0.990705621379;

	leb_tmp[138].p_cart.x = -0.096183085230;
	leb_tmp[138].p_cart.y = -0.096183085230;
	leb_tmp[138].p_cart.z = 0.990705621379;

	leb_tmp[139].p_cart.x = -0.096183085230;
	leb_tmp[139].p_cart.y = 0.096183085230;
	leb_tmp[139].p_cart.z = -0.990705621379;

	leb_tmp[140].p_cart.x = 0.096183085230;
	leb_tmp[140].p_cart.y = -0.096183085230;
	leb_tmp[140].p_cart.z = -0.990705621379;

	leb_tmp[141].p_cart.x = -0.096183085230;
	leb_tmp[141].p_cart.y = -0.096183085230;
	leb_tmp[141].p_cart.z = -0.990705621379;

	leb_tmp[142].p_cart.x = 0.096183085230;
	leb_tmp[142].p_cart.y = 0.990705621379;
	leb_tmp[142].p_cart.z = 0.096183085230;

	leb_tmp[143].p_cart.x = -0.096183085230;
	leb_tmp[143].p_cart.y = 0.990705621379;
	leb_tmp[143].p_cart.z = 0.096183085230;

	leb_tmp[144].p_cart.x = 0.096183085230;
	leb_tmp[144].p_cart.y = -0.990705621379;
	leb_tmp[144].p_cart.z = 0.096183085230;

	leb_tmp[145].p_cart.x = 0.096183085230;
	leb_tmp[145].p_cart.y = 0.990705621379;
	leb_tmp[145].p_cart.z = -0.096183085230;

	leb_tmp[146].p_cart.x = -0.096183085230;
	leb_tmp[146].p_cart.y = -0.990705621379;
	leb_tmp[146].p_cart.z = 0.096183085230;

	leb_tmp[147].p_cart.x = -0.096183085230;
	leb_tmp[147].p_cart.y = 0.990705621379;
	leb_tmp[147].p_cart.z = -0.096183085230;

	leb_tmp[148].p_cart.x = 0.096183085230;
	leb_tmp[148].p_cart.y = -0.990705621379;
	leb_tmp[148].p_cart.z = -0.096183085230;

	leb_tmp[149].p_cart.x = -0.096183085230;
	leb_tmp[149].p_cart.y = -0.990705621379;
	leb_tmp[149].p_cart.z = -0.096183085230;

	leb_tmp[150].p_cart.x = 0.990705621379;
	leb_tmp[150].p_cart.y = 0.096183085230;
	leb_tmp[150].p_cart.z = 0.096183085230;

	leb_tmp[151].p_cart.x = -0.990705621379;
	leb_tmp[151].p_cart.y = 0.096183085230;
	leb_tmp[151].p_cart.z = 0.096183085230;

	leb_tmp[152].p_cart.x = 0.990705621379;
	leb_tmp[152].p_cart.y = -0.096183085230;
	leb_tmp[152].p_cart.z = 0.096183085230;

	leb_tmp[153].p_cart.x = 0.990705621379;
	leb_tmp[153].p_cart.y = 0.096183085230;
	leb_tmp[153].p_cart.z = -0.096183085230;

	leb_tmp[154].p_cart.x = -0.990705621379;
	leb_tmp[154].p_cart.y = -0.096183085230;
	leb_tmp[154].p_cart.z = 0.096183085230;

	leb_tmp[155].p_cart.x = -0.990705621379;
	leb_tmp[155].p_cart.y = 0.096183085230;
	leb_tmp[155].p_cart.z = -0.096183085230;

	leb_tmp[156].p_cart.x = 0.990705621379;
	leb_tmp[156].p_cart.y = -0.096183085230;
	leb_tmp[156].p_cart.z = -0.096183085230;

	leb_tmp[157].p_cart.x = -0.990705621379;
	leb_tmp[157].p_cart.y = -0.096183085230;
	leb_tmp[157].p_cart.z = -0.096183085230;

	leb_tmp[158].p_cart.x = 0.820326419828;
	leb_tmp[158].p_cart.y = 0.571895589188;
	leb_tmp[158].p_cart.z = 0.000000000000;

	leb_tmp[159].p_cart.x = -0.820326419828;
	leb_tmp[159].p_cart.y = 0.571895589188;
	leb_tmp[159].p_cart.z = 0.000000000000;

	leb_tmp[160].p_cart.x = 0.820326419828;
	leb_tmp[160].p_cart.y = -0.571895589188;
	leb_tmp[160].p_cart.z = 0.000000000000;

	leb_tmp[161].p_cart.x = -0.820326419828;
	leb_tmp[161].p_cart.y = -0.571895589188;
	leb_tmp[161].p_cart.z = 0.000000000000;

	leb_tmp[162].p_cart.x = 0.571895589188;
	leb_tmp[162].p_cart.y = 0.820326419828;
	leb_tmp[162].p_cart.z = 0.000000000000;

	leb_tmp[163].p_cart.x = -0.571895589188;
	leb_tmp[163].p_cart.y = 0.820326419828;
	leb_tmp[163].p_cart.z = 0.000000000000;

	leb_tmp[164].p_cart.x = 0.571895589188;
	leb_tmp[164].p_cart.y = -0.820326419828;
	leb_tmp[164].p_cart.z = 0.000000000000;

	leb_tmp[165].p_cart.x = -0.571895589188;
	leb_tmp[165].p_cart.y = -0.820326419828;
	leb_tmp[165].p_cart.z = 0.000000000000;

	leb_tmp[166].p_cart.x = 0.820326419828;
	leb_tmp[166].p_cart.y = 0.000000000000;
	leb_tmp[166].p_cart.z = 0.571895589188;

	leb_tmp[167].p_cart.x = -0.820326419828;
	leb_tmp[167].p_cart.y = 0.000000000000;
	leb_tmp[167].p_cart.z = 0.571895589188;

	leb_tmp[168].p_cart.x = 0.820326419828;
	leb_tmp[168].p_cart.y = 0.000000000000;
	leb_tmp[168].p_cart.z = -0.571895589188;

	leb_tmp[169].p_cart.x = -0.820326419828;
	leb_tmp[169].p_cart.y = 0.000000000000;
	leb_tmp[169].p_cart.z = -0.571895589188;

	leb_tmp[170].p_cart.x = 0.571895589188;
	leb_tmp[170].p_cart.y = 0.000000000000;
	leb_tmp[170].p_cart.z = 0.820326419828;

	leb_tmp[171].p_cart.x = -0.571895589188;
	leb_tmp[171].p_cart.y = 0.000000000000;
	leb_tmp[171].p_cart.z = 0.820326419828;

	leb_tmp[172].p_cart.x = 0.571895589188;
	leb_tmp[172].p_cart.y = 0.000000000000;
	leb_tmp[172].p_cart.z = -0.820326419828;

	leb_tmp[173].p_cart.x = -0.571895589188;
	leb_tmp[173].p_cart.y = 0.000000000000;
	leb_tmp[173].p_cart.z = -0.820326419828;

	leb_tmp[174].p_cart.x = 0.000000000000;
	leb_tmp[174].p_cart.y = 0.820326419828;
	leb_tmp[174].p_cart.z = 0.571895589188;

	leb_tmp[175].p_cart.x = 0.000000000000;
	leb_tmp[175].p_cart.y = -0.820326419828;
	leb_tmp[175].p_cart.z = 0.571895589188;

	leb_tmp[176].p_cart.x = 0.000000000000;
	leb_tmp[176].p_cart.y = 0.820326419828;
	leb_tmp[176].p_cart.z = -0.571895589188;

	leb_tmp[177].p_cart.x = 0.000000000000;
	leb_tmp[177].p_cart.y = -0.820326419828;
	leb_tmp[177].p_cart.z = -0.571895589188;

	leb_tmp[178].p_cart.x = 0.000000000000;
	leb_tmp[178].p_cart.y = 0.571895589188;
	leb_tmp[178].p_cart.z = 0.820326419828;

	leb_tmp[179].p_cart.x = 0.000000000000;
	leb_tmp[179].p_cart.y = -0.571895589188;
	leb_tmp[179].p_cart.z = 0.820326419828;

	leb_tmp[180].p_cart.x = 0.000000000000;
	leb_tmp[180].p_cart.y = 0.571895589188;
	leb_tmp[180].p_cart.z = -0.820326419828;

	leb_tmp[181].p_cart.x = 0.000000000000;
	leb_tmp[181].p_cart.y = -0.571895589188;
	leb_tmp[181].p_cart.z = -0.820326419828;

	leb_tmp[182].p_cart.x = 0.964408914879;
	leb_tmp[182].p_cart.y = 0.264415288706;
	leb_tmp[182].p_cart.z = 0.000000000000;

	leb_tmp[183].p_cart.x = -0.964408914879;
	leb_tmp[183].p_cart.y = 0.264415288706;
	leb_tmp[183].p_cart.z = 0.000000000000;

	leb_tmp[184].p_cart.x = 0.964408914879;
	leb_tmp[184].p_cart.y = -0.264415288706;
	leb_tmp[184].p_cart.z = 0.000000000000;

	leb_tmp[185].p_cart.x = -0.964408914879;
	leb_tmp[185].p_cart.y = -0.264415288706;
	leb_tmp[185].p_cart.z = 0.000000000000;

	leb_tmp[186].p_cart.x = 0.264415288706;
	leb_tmp[186].p_cart.y = 0.964408914879;
	leb_tmp[186].p_cart.z = 0.000000000000;

	leb_tmp[187].p_cart.x = -0.264415288706;
	leb_tmp[187].p_cart.y = 0.964408914879;
	leb_tmp[187].p_cart.z = 0.000000000000;

	leb_tmp[188].p_cart.x = 0.264415288706;
	leb_tmp[188].p_cart.y = -0.964408914879;
	leb_tmp[188].p_cart.z = 0.000000000000;

	leb_tmp[189].p_cart.x = -0.264415288706;
	leb_tmp[189].p_cart.y = -0.964408914879;
	leb_tmp[189].p_cart.z = 0.000000000000;

	leb_tmp[190].p_cart.x = 0.964408914879;
	leb_tmp[190].p_cart.y = 0.000000000000;
	leb_tmp[190].p_cart.z = 0.264415288706;

	leb_tmp[191].p_cart.x = -0.964408914879;
	leb_tmp[191].p_cart.y = 0.000000000000;
	leb_tmp[191].p_cart.z = 0.264415288706;

	leb_tmp[192].p_cart.x = 0.964408914879;
	leb_tmp[192].p_cart.y = 0.000000000000;
	leb_tmp[192].p_cart.z = -0.264415288706;

	leb_tmp[193].p_cart.x = -0.964408914879;
	leb_tmp[193].p_cart.y = 0.000000000000;
	leb_tmp[193].p_cart.z = -0.264415288706;

	leb_tmp[194].p_cart.x = 0.264415288706;
	leb_tmp[194].p_cart.y = 0.000000000000;
	leb_tmp[194].p_cart.z = 0.964408914879;

	leb_tmp[195].p_cart.x = -0.264415288706;
	leb_tmp[195].p_cart.y = 0.000000000000;
	leb_tmp[195].p_cart.z = 0.964408914879;

	leb_tmp[196].p_cart.x = 0.264415288706;
	leb_tmp[196].p_cart.y = 0.000000000000;
	leb_tmp[196].p_cart.z = -0.964408914879;

	leb_tmp[197].p_cart.x = -0.264415288706;
	leb_tmp[197].p_cart.y = 0.000000000000;
	leb_tmp[197].p_cart.z = -0.964408914879;

	leb_tmp[198].p_cart.x = 0.000000000000;
	leb_tmp[198].p_cart.y = 0.964408914879;
	leb_tmp[198].p_cart.z = 0.264415288706;

	leb_tmp[199].p_cart.x = 0.000000000000;
	leb_tmp[199].p_cart.y = -0.964408914879;
	leb_tmp[199].p_cart.z = 0.264415288706;

	leb_tmp[200].p_cart.x = 0.000000000000;
	leb_tmp[200].p_cart.y = 0.964408914879;
	leb_tmp[200].p_cart.z = -0.264415288706;

	leb_tmp[201].p_cart.x = 0.000000000000;
	leb_tmp[201].p_cart.y = -0.964408914879;
	leb_tmp[201].p_cart.z = -0.264415288706;

	leb_tmp[202].p_cart.x = 0.000000000000;
	leb_tmp[202].p_cart.y = 0.264415288706;
	leb_tmp[202].p_cart.z = 0.964408914879;

	leb_tmp[203].p_cart.x = 0.000000000000;
	leb_tmp[203].p_cart.y = -0.264415288706;
	leb_tmp[203].p_cart.z = 0.964408914879;

	leb_tmp[204].p_cart.x = 0.000000000000;
	leb_tmp[204].p_cart.y = 0.264415288706;
	leb_tmp[204].p_cart.z = -0.964408914879;

	leb_tmp[205].p_cart.x = 0.000000000000;
	leb_tmp[205].p_cart.y = -0.264415288706;
	leb_tmp[205].p_cart.z = -0.964408914879;

	leb_tmp[206].p_cart.x = 0.251003475177;
	leb_tmp[206].p_cart.y = 0.800072749407;
	leb_tmp[206].p_cart.z = 0.544867737258;

	leb_tmp[207].p_cart.x = -0.251003475177;
	leb_tmp[207].p_cart.y = 0.800072749407;
	leb_tmp[207].p_cart.z = 0.544867737258;

	leb_tmp[208].p_cart.x = 0.251003475177;
	leb_tmp[208].p_cart.y = -0.800072749407;
	leb_tmp[208].p_cart.z = 0.544867737258;

	leb_tmp[209].p_cart.x = 0.251003475177;
	leb_tmp[209].p_cart.y = 0.800072749407;
	leb_tmp[209].p_cart.z = -0.544867737258;

	leb_tmp[210].p_cart.x = -0.251003475177;
	leb_tmp[210].p_cart.y = -0.800072749407;
	leb_tmp[210].p_cart.z = 0.544867737258;

	leb_tmp[211].p_cart.x = -0.251003475177;
	leb_tmp[211].p_cart.y = 0.800072749407;
	leb_tmp[211].p_cart.z = -0.544867737258;

	leb_tmp[212].p_cart.x = 0.251003475177;
	leb_tmp[212].p_cart.y = -0.800072749407;
	leb_tmp[212].p_cart.z = -0.544867737258;

	leb_tmp[213].p_cart.x = -0.251003475177;
	leb_tmp[213].p_cart.y = -0.800072749407;
	leb_tmp[213].p_cart.z = -0.544867737258;

	leb_tmp[214].p_cart.x = 0.251003475177;
	leb_tmp[214].p_cart.y = 0.544867737258;
	leb_tmp[214].p_cart.z = 0.800072749407;

	leb_tmp[215].p_cart.x = -0.251003475177;
	leb_tmp[215].p_cart.y = 0.544867737258;
	leb_tmp[215].p_cart.z = 0.800072749407;

	leb_tmp[216].p_cart.x = 0.251003475177;
	leb_tmp[216].p_cart.y = -0.544867737258;
	leb_tmp[216].p_cart.z = 0.800072749407;

	leb_tmp[217].p_cart.x = 0.251003475177;
	leb_tmp[217].p_cart.y = 0.544867737258;
	leb_tmp[217].p_cart.z = -0.800072749407;

	leb_tmp[218].p_cart.x = -0.251003475177;
	leb_tmp[218].p_cart.y = -0.544867737258;
	leb_tmp[218].p_cart.z = 0.800072749407;

	leb_tmp[219].p_cart.x = -0.251003475177;
	leb_tmp[219].p_cart.y = 0.544867737258;
	leb_tmp[219].p_cart.z = -0.800072749407;

	leb_tmp[220].p_cart.x = 0.251003475177;
	leb_tmp[220].p_cart.y = -0.544867737258;
	leb_tmp[220].p_cart.z = -0.800072749407;

	leb_tmp[221].p_cart.x = -0.251003475177;
	leb_tmp[221].p_cart.y = -0.544867737258;
	leb_tmp[221].p_cart.z = -0.800072749407;

	leb_tmp[222].p_cart.x = 0.800072749407;
	leb_tmp[222].p_cart.y = 0.251003475177;
	leb_tmp[222].p_cart.z = 0.544867737258;

	leb_tmp[223].p_cart.x = -0.800072749407;
	leb_tmp[223].p_cart.y = 0.251003475177;
	leb_tmp[223].p_cart.z = 0.544867737258;

	leb_tmp[224].p_cart.x = 0.800072749407;
	leb_tmp[224].p_cart.y = -0.251003475177;
	leb_tmp[224].p_cart.z = 0.544867737258;

	leb_tmp[225].p_cart.x = 0.800072749407;
	leb_tmp[225].p_cart.y = 0.251003475177;
	leb_tmp[225].p_cart.z = -0.544867737258;

	leb_tmp[226].p_cart.x = -0.800072749407;
	leb_tmp[226].p_cart.y = -0.251003475177;
	leb_tmp[226].p_cart.z = 0.544867737258;

	leb_tmp[227].p_cart.x = -0.800072749407;
	leb_tmp[227].p_cart.y = 0.251003475177;
	leb_tmp[227].p_cart.z = -0.544867737258;

	leb_tmp[228].p_cart.x = 0.800072749407;
	leb_tmp[228].p_cart.y = -0.251003475177;
	leb_tmp[228].p_cart.z = -0.544867737258;

	leb_tmp[229].p_cart.x = -0.800072749407;
	leb_tmp[229].p_cart.y = -0.251003475177;
	leb_tmp[229].p_cart.z = -0.544867737258;

	leb_tmp[230].p_cart.x = 0.800072749407;
	leb_tmp[230].p_cart.y = 0.544867737258;
	leb_tmp[230].p_cart.z = 0.251003475177;

	leb_tmp[231].p_cart.x = -0.800072749407;
	leb_tmp[231].p_cart.y = 0.544867737258;
	leb_tmp[231].p_cart.z = 0.251003475177;

	leb_tmp[232].p_cart.x = 0.800072749407;
	leb_tmp[232].p_cart.y = -0.544867737258;
	leb_tmp[232].p_cart.z = 0.251003475177;

	leb_tmp[233].p_cart.x = 0.800072749407;
	leb_tmp[233].p_cart.y = 0.544867737258;
	leb_tmp[233].p_cart.z = -0.251003475177;

	leb_tmp[234].p_cart.x = -0.800072749407;
	leb_tmp[234].p_cart.y = -0.544867737258;
	leb_tmp[234].p_cart.z = 0.251003475177;

	leb_tmp[235].p_cart.x = -0.800072749407;
	leb_tmp[235].p_cart.y = 0.544867737258;
	leb_tmp[235].p_cart.z = -0.251003475177;

	leb_tmp[236].p_cart.x = 0.800072749407;
	leb_tmp[236].p_cart.y = -0.544867737258;
	leb_tmp[236].p_cart.z = -0.251003475177;

	leb_tmp[237].p_cart.x = -0.800072749407;
	leb_tmp[237].p_cart.y = -0.544867737258;
	leb_tmp[237].p_cart.z = -0.251003475177;

	leb_tmp[238].p_cart.x = 0.544867737258;
	leb_tmp[238].p_cart.y = 0.251003475177;
	leb_tmp[238].p_cart.z = 0.800072749407;

	leb_tmp[239].p_cart.x = -0.544867737258;
	leb_tmp[239].p_cart.y = 0.251003475177;
	leb_tmp[239].p_cart.z = 0.800072749407;

	leb_tmp[240].p_cart.x = 0.544867737258;
	leb_tmp[240].p_cart.y = -0.251003475177;
	leb_tmp[240].p_cart.z = 0.800072749407;

	leb_tmp[241].p_cart.x = 0.544867737258;
	leb_tmp[241].p_cart.y = 0.251003475177;
	leb_tmp[241].p_cart.z = -0.800072749407;

	leb_tmp[242].p_cart.x = -0.544867737258;
	leb_tmp[242].p_cart.y = -0.251003475177;
	leb_tmp[242].p_cart.z = 0.800072749407;

	leb_tmp[243].p_cart.x = -0.544867737258;
	leb_tmp[243].p_cart.y = 0.251003475177;
	leb_tmp[243].p_cart.z = -0.800072749407;

	leb_tmp[244].p_cart.x = 0.544867737258;
	leb_tmp[244].p_cart.y = -0.251003475177;
	leb_tmp[244].p_cart.z = -0.800072749407;

	leb_tmp[245].p_cart.x = -0.544867737258;
	leb_tmp[245].p_cart.y = -0.251003475177;
	leb_tmp[245].p_cart.z = -0.800072749407;

	leb_tmp[246].p_cart.x = 0.544867737258;
	leb_tmp[246].p_cart.y = 0.800072749407;
	leb_tmp[246].p_cart.z = 0.251003475177;

	leb_tmp[247].p_cart.x = -0.544867737258;
	leb_tmp[247].p_cart.y = 0.800072749407;
	leb_tmp[247].p_cart.z = 0.251003475177;

	leb_tmp[248].p_cart.x = 0.544867737258;
	leb_tmp[248].p_cart.y = -0.800072749407;
	leb_tmp[248].p_cart.z = 0.251003475177;

	leb_tmp[249].p_cart.x = 0.544867737258;
	leb_tmp[249].p_cart.y = 0.800072749407;
	leb_tmp[249].p_cart.z = -0.251003475177;

	leb_tmp[250].p_cart.x = -0.544867737258;
	leb_tmp[250].p_cart.y = -0.800072749407;
	leb_tmp[250].p_cart.z = 0.251003475177;

	leb_tmp[251].p_cart.x = -0.544867737258;
	leb_tmp[251].p_cart.y = 0.800072749407;
	leb_tmp[251].p_cart.z = -0.251003475177;

	leb_tmp[252].p_cart.x = 0.544867737258;
	leb_tmp[252].p_cart.y = -0.800072749407;
	leb_tmp[252].p_cart.z = -0.251003475177;

	leb_tmp[253].p_cart.x = -0.544867737258;
	leb_tmp[253].p_cart.y = -0.800072749407;
	leb_tmp[253].p_cart.z = -0.251003475177;

	leb_tmp[254].p_cart.x = 0.902442529533;
	leb_tmp[254].p_cart.y = 0.412772408317;
	leb_tmp[254].p_cart.z = 0.123354853258;

	leb_tmp[255].p_cart.x = -0.902442529533;
	leb_tmp[255].p_cart.y = 0.412772408317;
	leb_tmp[255].p_cart.z = 0.123354853258;

	leb_tmp[256].p_cart.x = 0.902442529533;
	leb_tmp[256].p_cart.y = -0.412772408317;
	leb_tmp[256].p_cart.z = 0.123354853258;

	leb_tmp[257].p_cart.x = 0.902442529533;
	leb_tmp[257].p_cart.y = 0.412772408317;
	leb_tmp[257].p_cart.z = -0.123354853258;

	leb_tmp[258].p_cart.x = -0.902442529533;
	leb_tmp[258].p_cart.y = -0.412772408317;
	leb_tmp[258].p_cart.z = 0.123354853258;

	leb_tmp[259].p_cart.x = -0.902442529533;
	leb_tmp[259].p_cart.y = 0.412772408317;
	leb_tmp[259].p_cart.z = -0.123354853258;

	leb_tmp[260].p_cart.x = 0.902442529533;
	leb_tmp[260].p_cart.y = -0.412772408317;
	leb_tmp[260].p_cart.z = -0.123354853258;

	leb_tmp[261].p_cart.x = -0.902442529533;
	leb_tmp[261].p_cart.y = -0.412772408317;
	leb_tmp[261].p_cart.z = -0.123354853258;

	leb_tmp[262].p_cart.x = 0.902442529533;
	leb_tmp[262].p_cart.y = 0.123354853258;
	leb_tmp[262].p_cart.z = 0.412772408317;

	leb_tmp[263].p_cart.x = -0.902442529533;
	leb_tmp[263].p_cart.y = 0.123354853258;
	leb_tmp[263].p_cart.z = 0.412772408317;

	leb_tmp[264].p_cart.x = 0.902442529533;
	leb_tmp[264].p_cart.y = -0.123354853258;
	leb_tmp[264].p_cart.z = 0.412772408317;

	leb_tmp[265].p_cart.x = 0.902442529533;
	leb_tmp[265].p_cart.y = 0.123354853258;
	leb_tmp[265].p_cart.z = -0.412772408317;

	leb_tmp[266].p_cart.x = -0.902442529533;
	leb_tmp[266].p_cart.y = -0.123354853258;
	leb_tmp[266].p_cart.z = 0.412772408317;

	leb_tmp[267].p_cart.x = -0.902442529533;
	leb_tmp[267].p_cart.y = 0.123354853258;
	leb_tmp[267].p_cart.z = -0.412772408317;

	leb_tmp[268].p_cart.x = 0.902442529533;
	leb_tmp[268].p_cart.y = -0.123354853258;
	leb_tmp[268].p_cart.z = -0.412772408317;

	leb_tmp[269].p_cart.x = -0.902442529533;
	leb_tmp[269].p_cart.y = -0.123354853258;
	leb_tmp[269].p_cart.z = -0.412772408317;

	leb_tmp[270].p_cart.x = 0.412772408317;
	leb_tmp[270].p_cart.y = 0.902442529533;
	leb_tmp[270].p_cart.z = 0.123354853258;

	leb_tmp[271].p_cart.x = -0.412772408317;
	leb_tmp[271].p_cart.y = 0.902442529533;
	leb_tmp[271].p_cart.z = 0.123354853258;

	leb_tmp[272].p_cart.x = 0.412772408317;
	leb_tmp[272].p_cart.y = -0.902442529533;
	leb_tmp[272].p_cart.z = 0.123354853258;

	leb_tmp[273].p_cart.x = 0.412772408317;
	leb_tmp[273].p_cart.y = 0.902442529533;
	leb_tmp[273].p_cart.z = -0.123354853258;

	leb_tmp[274].p_cart.x = -0.412772408317;
	leb_tmp[274].p_cart.y = -0.902442529533;
	leb_tmp[274].p_cart.z = 0.123354853258;

	leb_tmp[275].p_cart.x = -0.412772408317;
	leb_tmp[275].p_cart.y = 0.902442529533;
	leb_tmp[275].p_cart.z = -0.123354853258;

	leb_tmp[276].p_cart.x = 0.412772408317;
	leb_tmp[276].p_cart.y = -0.902442529533;
	leb_tmp[276].p_cart.z = -0.123354853258;

	leb_tmp[277].p_cart.x = -0.412772408317;
	leb_tmp[277].p_cart.y = -0.902442529533;
	leb_tmp[277].p_cart.z = -0.123354853258;

	leb_tmp[278].p_cart.x = 0.412772408317;
	leb_tmp[278].p_cart.y = 0.123354853258;
	leb_tmp[278].p_cart.z = 0.902442529533;

	leb_tmp[279].p_cart.x = -0.412772408317;
	leb_tmp[279].p_cart.y = 0.123354853258;
	leb_tmp[279].p_cart.z = 0.902442529533;

	leb_tmp[280].p_cart.x = 0.412772408317;
	leb_tmp[280].p_cart.y = -0.123354853258;
	leb_tmp[280].p_cart.z = 0.902442529533;
	
	leb_tmp[281].p_cart.x = 0.412772408317;
	leb_tmp[281].p_cart.y = 0.123354853258;
	leb_tmp[281].p_cart.z = -0.902442529533;

	leb_tmp[282].p_cart.x = -0.412772408317;
	leb_tmp[282].p_cart.y = -0.123354853258;
	leb_tmp[282].p_cart.z = 0.902442529533;

	leb_tmp[283].p_cart.x = -0.412772408317;
	leb_tmp[283].p_cart.y = 0.123354853258;
	leb_tmp[283].p_cart.z = -0.902442529533;

	leb_tmp[284].p_cart.x = 0.412772408317;
	leb_tmp[284].p_cart.y = -0.123354853258;
	leb_tmp[284].p_cart.z = -0.902442529533;

	leb_tmp[285].p_cart.x = -0.412772408317;
	leb_tmp[285].p_cart.y = -0.123354853258;
	leb_tmp[285].p_cart.z = -0.902442529533;

	leb_tmp[286].p_cart.x = 0.123354853258;
	leb_tmp[286].p_cart.y = 0.902442529533;
	leb_tmp[286].p_cart.z = 0.412772408317;

	leb_tmp[287].p_cart.x = -0.123354853258;
	leb_tmp[287].p_cart.y = 0.902442529533;
	leb_tmp[287].p_cart.z = 0.412772408317;

	leb_tmp[288].p_cart.x = 0.123354853258;
	leb_tmp[288].p_cart.y = -0.902442529533;
	leb_tmp[288].p_cart.z = 0.412772408317;

	leb_tmp[289].p_cart.x = 0.123354853258;
	leb_tmp[289].p_cart.y = 0.902442529533;
	leb_tmp[289].p_cart.z = -0.412772408317;

	leb_tmp[290].p_cart.x = -0.123354853258;
	leb_tmp[290].p_cart.y = -0.902442529533;
	leb_tmp[290].p_cart.z = 0.412772408317;

	leb_tmp[291].p_cart.x = -0.123354853258;
	leb_tmp[291].p_cart.y = 0.902442529533;
	leb_tmp[291].p_cart.z = -0.412772408317;

	leb_tmp[292].p_cart.x = 0.123354853258;
	leb_tmp[292].p_cart.y = -0.902442529533;
	leb_tmp[292].p_cart.z = -0.412772408317;

	leb_tmp[293].p_cart.x = -0.123354853258;
	leb_tmp[293].p_cart.y = -0.902442529533;
	leb_tmp[293].p_cart.z = -0.412772408317;

	leb_tmp[294].p_cart.x = 0.123354853258;
	leb_tmp[294].p_cart.y = 0.412772408317;
	leb_tmp[294].p_cart.z = 0.902442529533;

	leb_tmp[295].p_cart.x = -0.123354853258;
	leb_tmp[295].p_cart.y = 0.412772408317;
	leb_tmp[295].p_cart.z = 0.902442529533;

	leb_tmp[296].p_cart.x = 0.123354853258;
	leb_tmp[296].p_cart.y = -0.412772408317;
	leb_tmp[296].p_cart.z = 0.902442529533;

	leb_tmp[297].p_cart.x = 0.123354853258;
	leb_tmp[297].p_cart.y = 0.412772408317;
	leb_tmp[297].p_cart.z = -0.902442529533;

	leb_tmp[298].p_cart.x = -0.123354853258;
	leb_tmp[298].p_cart.y = -0.412772408317;
	leb_tmp[298].p_cart.z = 0.902442529533;

	leb_tmp[299].p_cart.x = -0.123354853258;
	leb_tmp[299].p_cart.y = 0.412772408317;
	leb_tmp[299].p_cart.z = -0.902442529533;

	leb_tmp[300].p_cart.x = 0.123354853258;
	leb_tmp[300].p_cart.y = -0.412772408317;
	leb_tmp[300].p_cart.z = -0.902442529533;

	leb_tmp[301].p_cart.x = -0.123354853258;
	leb_tmp[301].p_cart.y = -0.412772408317;
	leb_tmp[301].p_cart.z = -0.902442529533;
	break;
    
    default:
	punt("Angular grid specified not implemented");
    }
    
    return leb_tmp;
}

