/*****************************************************************************

	INTERNALS.CC

	member functions for internals class
******************************************************************************/


extern "C" {
  #include <stdio.h>
  #include <stdlib.h>
  #include <string.h>
  #include <math.h>
  #include <libciomr.h>
  #include <ip_libv1.h>
  #include <physconst.h>
  #include "bond_lengths.h"
}

#define EXTERN
#include "opt.h"
#include "cartesians.h"
#include "internals.h"

// make a<b in stretch a-b, a<c in bend a-b-c and a<d in tors a-b-c-d
extern void swap_tors(int *a, int *b, int *c, int *d);
extern void swap(int *a, int *b);



/*-----------------------------------------------------------------------------

	INTERNALS::INDEX_TO_ID

	Converts the absolute optking index into the user-assigned id number 
-----------------------------------------------------------------------------*/

int internals :: index_to_id(int index) {
   int i, count = -1, id = -1;
   for (i=0;i<stre.get_num();++i)
      if (++count == index) {id = stre.get_id(i); break;}
   for (i=0;i<bend.get_num();++i)
      if (++count == index) {id = bend.get_id(i); break;}
   for (i=0;i<tors.get_num();++i)
      if (++count == index) {id = tors.get_id(i); break;}
   for (i=0;i<out.get_num();++i)
      if (++count == index) {id = out.get_id(i); break;}
//   fprintf(outfile,"index_to_id(%d): returning id = %d\n",index, id);
   return id;
}



/*-----------------------------------------------------------------------------

	INTERNALS::ID_TO_INDEX

	Converts the user-assigned id number into an absolute optking index 
-----------------------------------------------------------------------------*/

int internals :: id_to_index(int id) {
   int i, count = 0, index = -1, found = 0;
   for (i=0;i<stre.get_num();++i,++count)
      if (id == stre.get_id(i)) {index = count; ++found;}
   for (i=0;i<bend.get_num();++i,++count)
      if (id == bend.get_id(i)) {index = count; ++found;}
   for (i=0;i<tors.get_num();++i,++count)
      if (id == tors.get_id(i)) {index = count; ++found;}
   for (i=0;i<out.get_num();++i,++count)
      if (id == out.get_id(i)) {index = count; ++found;}
   if (found == 0) {
      fprintf(outfile,"Error: Simple internal id not found.\n");
      exit(2);
   }
   else if (found > 1) {
      fprintf(outfile,"Error: Simple internal id appears more than once.\n");
      exit(2);
   }
// fprintf(outfile,"id_to_index(%d): returning index = %d\n",id, index);
   return index;
}



/*----------------------------------------------------------------------------

	INTERNALS::LOCATE_ID

	Given the user-assigned id number, returns the type of internal
	and the subindex within that type 
----------------------------------------------------------------------------*/

void internals :: locate_id(int id, int *intco_type, int *sub_index) {
   int i, found = 0;
// fprintf(outfile,"locating id %d\n",id);
   for (i=0;i<stre.get_num();++i)
      if (stre.get_id(i) == id) {
         ++found; *intco_type = STRE_TYPE; *sub_index = i;
      }
   for (i=0;i<bend.get_num();++i)
      if (bend.get_id(i) == id) {
         ++found; *intco_type = BEND_TYPE; *sub_index = i;
      }
   for (i=0;i<tors.get_num();++i)
      if (tors.get_id(i) == id) {
         ++found; *intco_type = TORS_TYPE; *sub_index = i;
      }
   for (i=0;i<out.get_num();++i)
      if (out.get_id(i) == id) {
         ++found; *intco_type = OUT_TYPE; *sub_index = i;
      }
   if (found == 0) {
     fprintf(outfile,"Error: Simple internal id not found.\n");
     exit(2);
   }
   else if (found > 1) {
     fprintf(outfile,"Error: Simple internal id appears more than once.\n");
     exit(2);
   }
// fprintf(outfile,"returning intco_type: %d sub_index %d\n",*intco_type, *sub_index);
   return;
}




/*??????????*/

internals :: ~internals() {
  fprintf(outfile,"destructor for internals called\n");
  fflush(outfile);
// Destructors are called automatically ??
//  stre.~stretch_set();
//  bend.~bend_set();
//  tors.~torsion_set();
//  out.~out_set();
}





/*-----------------------------------------------------------------------------

	INTERNALS::INTERNALS

	class constructor for internals
-----------------------------------------------------------------------------*/

internals :: internals(cartesians& carts, int user_intcos, int *size_arr)
  : stre(size_arr[0]), bend(size_arr[1]), tors(size_arr[2]), out(size_arr[3]) 
{

 if (user_intcos) {
  /* Read in simple internal coordinates from intco.dat */
  int i,j,a,b,c,d,e;
  fp_intco = fopen("intco.dat","r");
  ip_set_uppercase(1);
  ip_initialize(fp_intco, outfile);
  ip_cwk_add(":INTCO");

  i=0;
  if (ip_exist("STRE",0)) {
     ip_count("STRE",&i,0);
     stre.set_num(i);

     for(i=0;i<stre.get_num();++i) {
        ip_count("STRE",&j,1,i);
        if (j != 3) {
           fprintf(outfile,"Stretch %d is of wrong dimension.\n",i+1);
           exit(2);
        }
        ip_data("STRE","%d",&(a),2,i,0);
        ip_data("STRE","%d",&(b),2,i,1);
        ip_data("STRE","%d",&(c),2,i,2);
        stre.set_id(i,a);
        swap(&b,&c);
        stre.set_A(i,b-1);
        stre.set_B(i,c-1);
        stre.set_val(i,0.0);
     }
  }
  else if (!ip_exist("STRETCH",0)) {stre.set_num(0);}
  
  i = 0;
  if (ip_exist("BEND",0)) {
     ip_count("BEND",&i,0);
     bend.set_num(i);

     for(i=0;i<bend.get_num();++i) {
        ip_count("BEND",&j,1,i);
        if (j != 4) {
           fprintf(outfile,"Bend %d is of wrong dimension.\n",i+1);
           exit(2);
        }
        ip_data("BEND","%d",&(a),2,i,0);  ip_data("BEND","%d",&(b),2,i,1);
        ip_data("BEND","%d",&(c),2,i,2);  ip_data("BEND","%d",&(d),2,i,3);
        bend.set_id(i,a);
        swap(&b,&d);
        bend.set_A(i,b-1);
        bend.set_B(i,c-1);
        bend.set_C(i,d-1);
        bend.set_val(i,0.0);
     }
  }
  else if (!ip_exist("BEND",0)) { bend.set_num(0);}
  
  i = 0;
  if (ip_exist("TORS",0)) {
     ip_count("TORS",&i,0);
     tors.set_num(i);

     for(i=0;i<tors.get_num();++i) {
        ip_count("TORS",&j,1,i);
        if (j != 5) {
           fprintf(outfile,"Torsion %d is of wrong dimension.\n",i+1);
           exit(2);
        }
        ip_data("TORS","%d",&(a),2,i,0);
        ip_data("TORS","%d",&(b),2,i,1);
        ip_data("TORS","%d",&(c),2,i,2);
        ip_data("TORS","%d",&(d),2,i,3);
        ip_data("TORS","%d",&(e),2,i,4);
        tors.set_id(i,a);
        swap_tors(&b,&c,&d,&e);
        tors.set_A(i,b-1);
        tors.set_B(i,c-1);
        tors.set_C(i,d-1);
        tors.set_D(i,e-1);
        tors.set_val(i,0.0);
     }
  }
  else if(!ip_exist("TORS",0)) { tors.set_num(0); }
  
  i = 0;
  if (ip_exist("OUT",0)) {
     ip_count("OUT",&i,0);
     out.set_num(i);

     for(i=0;i<out.get_num();++i) {
        ip_count("OUT",&j,1,i);
        if (j != 5) {
           fprintf(outfile,"Out-of-plane %d is of wrong dimension.\n",i+1);
           exit(2);
        }
        ip_data("OUT","%d",&(a),2,i,0);
        ip_data("OUT","%d",&(b),2,i,1);
        ip_data("OUT","%d",&(c),2,i,2);
        ip_data("OUT","%d",&(d),2,i,3);
        ip_data("OUT","%d",&(e),2,i,4);
        out.set_id(i,a);
        out.set_A(i,b-1);
        out.set_B(i,c-1);
        out.set_C(i,d-1);
        out.set_D(i,e-1);
        out.set_val(i,0.0);
      }
   }
  else if(!ip_exist("OUT",0)) out.set_num(0);
  
   ip_done();
 }
 else {

   /* Generate simple internal coordinates from a cartesian geometry */
  int i,j,k,a,b,c,d,e, *ioff, count, Zmax, Zmin, num_atoms;
  double *atom_dist, *coord;
  int **bonds, type_count=0, id_count=0;

  num_atoms = carts.get_num_atoms();
  coord = carts.get_coord();

  ioff = (int *) malloc (32641 * sizeof(int)) ;
  ioff[0] = 0 ;
  for (i = 1; i < 32641 ; i++)
    ioff[i] = ioff[i-1] + i;
 
  /* Compute atomic distance matrix */
  atom_dist = init_array( ((num_atoms+1)*num_atoms)/2 );
  count = -1;
  for (i=0;i<num_atoms;++i)
    for (j=0;j<=i;++j)
      atom_dist[++count]= sqrt(SQR(coord[3*i+0] - coord[3*j+0])+
                               SQR(coord[3*i+1] - coord[3*j+1])+
                               SQR(coord[3*i+2] - coord[3*j+2]));

  /* Determine which bonds are present */
  bonds = init_int_matrix(num_atoms, num_atoms);
  int count_of_bonds = 0;
  for (i=0;i<num_atoms;++i)
    for (j=0;j<i;++j) {
      Zmax = MAX((int)carts.get_atomic_num(i),(int)carts.get_atomic_num(j));
      Zmin = MIN((int)carts.get_atomic_num(i),(int)carts.get_atomic_num(j));
      a = ioff[Zmax-1] + (Zmin-1);
      if (bondl[a] != 0.0) {
        if (atom_dist[ioff[i]+j] < (1.2 * bondl[a])) {
           bonds[i][j] = 1;
           bonds[j][i] = 1;
        }
      }
      else  {
        fprintf(outfile,"WARNING! Optking does not know what bond lengths");
        fprintf(outfile,"to expect for all the atoms.\n");
        fprintf(outfile,"You may have to specify connectivity in input.");
      }
    }

  int num_bonds, num_nobonds;
  rewind(fp_input);
  ip_set_uppercase(1);
  ip_initialize(fp_input,outfile);
// ip_cwk_clear();
  ip_cwk_add(":OPTKING");
  if (!ip_exist("BONDS",0))
     num_bonds = 0;
  else
    ip_count("BONDS",&num_bonds,0,0);
  if (optinfo.print_simples)
    fprintf(outfile,"\nNumber of user-specified bonds: %3d\n",num_bonds);
  for(i=0;i<num_bonds;++i) {
    ip_data("BONDS","%d",&a,2,i,0);
    ip_data("BONDS","%d",&b,2,i,1);
    fprintf(outfile,"Bond %3d: %3d %3d\n",i+1,a,b);
    a -= 1;
    b -= 1;
    bonds[a][b] = 1;
    bonds[b][a] = 1;
  }
  if (!ip_exist("NOBONDS",0))
    num_nobonds = 0;
  else
    ip_count("NOBONDS",&num_nobonds,0,0);
  if (optinfo.print_simples)
    fprintf(outfile,"\nNumber of user-specified \"nobonds\": %3d\n",num_nobonds);
  for(i=0;i<num_nobonds;++i) {
     ip_data("NOBONDS","%d",&a,2,i,0);
     ip_data("NOBONDS","%d",&b,2,i,1);
     fprintf(outfile,"Nobond %3d: %3d %3d\n",i+1,a,b);
     a -= 1;
     b -= 1;
     bonds[a][b] = 0;
     bonds[b][a] = 0;
  }
  ip_done();

  fprintf(outfile,"\ngenerating internals\n"); fflush(outfile);
  id_count = 0;
  type_count = 0;
  for (i=0;i<num_atoms;++i)
    for (j=i+1;j<num_atoms;++j)
      if (bonds[i][j]) {
        stre.set_id(type_count,++id_count); // IDs set starting with 1
        stre.set_A(type_count,i);
        stre.set_B(type_count,j);
        stre.set_val(type_count++,0.0); // values computed later
      }
  stre.set_num(type_count);

  if (optinfo.print_simples) {
    fprintf(outfile,"\n+++ Bond Connectivity +++\n");
    for (i=0;i<num_atoms;++i) {
      fprintf(outfile,"%d:",i+1);
      for (j=0;j<num_atoms;++j)
        if (bonds[i][j]) fprintf(outfile," %d",j+1);
      fprintf(outfile,"\n");
    }
  }

  /* Determine bond angles */
  type_count = 0;
  for(i=0;i<num_atoms;++i)
    for(j=0;j<num_atoms;++j)
      for(k=i+1;k<num_atoms;++k)
        if (bonds[i][j] && bonds[j][k]) {
          bend.set_id(type_count,++id_count);
          bend.set_A(type_count,i);
          bend.set_B(type_count,j);
          bend.set_C(type_count,k);
          bend.set_val(type_count++,0.0);
        }
  bend.set_num(type_count);

  /* Determine torsions */
  type_count = 0;
  for(a=0;a<num_atoms;++a)
    for(b=0;b<num_atoms;++b)
      for(c=0;c<num_atoms;++c)
        if (c != a) {
          for(d=a+1;d<num_atoms;++d)
            if ( (d != b) && bonds[a][b] && bonds[b][c] && bonds[c][d]) {
              tors.set_id(type_count,++id_count);
              tors.set_A(type_count,a);
              tors.set_B(type_count,b);
              tors.set_C(type_count,c);
              tors.set_D(type_count,d);
              tors.set_val(type_count++,0.0);
            }
        }

   tors.set_num(type_count);
   out.set_num(0);
   fp_intco = fopen("intco.dat","w");
   print(fp_intco,0); 
 }

 fclose(fp_intco);
 num = stre.get_num() + bend.get_num() + tors.get_num() + out.get_num();
 if (num == 0) {
   fprintf(outfile,"Error: No simple internals were read.\n");
   exit(2);
 }

 return;

}





/*-----------------------------------------------------------------------------

	INTERNALS::PRINT

        Print out internal coordinates to a file in intco.dat format
	(print_flag == 0) or with values (print_flag == 1) 
----------------------------------------------------------------------------*/

/* Print out internal coordinates to a file in
   intco.dat format (print_flag == 0) or
   with values (print_flag == 1) */
void internals :: print(FILE *fp_out, int print_flag) {
  if (print_flag == 0) fprintf(fp_out, "intco: (\n");
  stre.print(fp_out, print_flag);
  bend.print(fp_out, print_flag);
  tors.print(fp_out, print_flag);
  out.print(fp_out, print_flag);
  if (print_flag == 0) fprintf(fp_out, ")\n");
  return;
}



/*----------------------------------------------------------------------------

	INTERNALS::COMPUTE
----------------------------------------------------------------------------*/

void internals :: compute_internals(int num_atoms, double *geom) {
  stre.compute(num_atoms, geom);
  bend.compute(num_atoms, geom);
  tors.compute(num_atoms, geom);
  out.compute(num_atoms, geom);
  return;
}



/*-----------------------------------------------------------------------------

	INTERNALS::COMPUTE_S
-----------------------------------------------------------------------------*/

void internals :: compute_s(int num_atoms, double *geom) {
  stre.compute_s(num_atoms, geom);
  bend.compute_s(num_atoms, geom);
  tors.compute_s(num_atoms, geom);
  out.compute_s(num_atoms, geom);
  return;
}



/*----------------------------------------------------------------------------

	INTERNALS::PRINT_S
----------------------------------------------------------------------------*/

void internals :: print_s() {
  stre.print_s();
  bend.print_s();
  tors.print_s();
  out.print_s();
  return;
}


