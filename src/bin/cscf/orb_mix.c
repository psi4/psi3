/*-----------------------------------------------------------

  orb_mix(); - A function to take the alpha and beta orbitals
  and destroy their alpha and beta symmetry to find unique 
  UHF solutions.

  The manner in which this will be done is to take the LUMO
  and make them mix them in the following way

  LUMO (alpha) = LUMO + HOMO
  LUMO (beta) = LUMO - HOMO

  Shawn Brown - (11/4/99)

  Revised on 15 April 2002 by David Sherrill
  This code made no sense before.  Orbitals were not
  kept orthonormal, uninitialized arrays were used, and
  numbering started from 1 instead of 0 (!).  Here is
  a less insane attempt.  
  ---------------------------------------------------------*/

#define EXTERN
#include "includes.h"
#include "common.h"
#ifndef PI
  #define PI 3.14159265359
#endif

void orb_mix(void)
{
    int i, num_mo, spin, irr, found_guess, homo_orb, lumo_orb, errcod;
    int lumo_irrep[2], homo_irrep[2];
    double lumo_energy[2], homo_energy[2];
    double c1, c2, tval, tval1, tval2, mixing_fraction;
 
    fprintf(outfile, "\n  ***** Orbital Mixing being used to find unique UHF solution *****\n");

    /* figure out which irrep has the HOMO for alpha spin */
    spin = 0;

    /* these guess eigenvalues are really awful but it seems ok
    for (irr=0; irr<num_ir; irr++) {
      fprintf(outfile, "Irrep %d\n", irr);
      for (i=0; i<scf_info[irr].num_mo; i++) {
        fprintf(outfile, "%lf\n", spin_info[spin].scf_spin[irr].fock_evals[i]);
      }
    }
    */

    for (irr=0,found_guess=0; irr<num_ir; irr++) {
      if (spin_info[spin].scf_spin[irr].noccup > 0) {
        tval = spin_info[spin].scf_spin[irr].fock_evals[spin_info[spin].scf_spin[irr].noccup-1]; 
        if (!found_guess) {
          homo_irrep[spin] = irr;
          homo_energy[spin] = tval; 
          found_guess = 1;
        }
        else {
          if (tval > homo_energy[spin]) {
            homo_energy[spin] = tval;
            homo_irrep[spin] = irr;
          }
        }
      } 
    } /* end loop over irreps */

    if (!found_guess) {
      fprintf(outfile, 
              "cscf: (orb_mix): Can't find a valid guess for HOMO irrep!\n");
      exit(1);
    } 


    /* now try to figure out where the LUMO is */

    for (irr=0,found_guess=0; irr<num_ir; irr++) {

      /* if there are any virtuals in this irrep... */
      if (scf_info[irr].num_mo > spin_info[spin].scf_spin[irr].noccup) {

        tval = spin_info[spin].scf_spin[irr].fock_evals[spin_info[spin].scf_spin[irr].noccup]; 
        if (!found_guess) {
          lumo_irrep[spin] = irr;
          lumo_energy[spin] = tval; 
          found_guess = 1;
        }
        else {
          if (tval < lumo_energy[spin]) {
            lumo_energy[spin] = tval;
            lumo_irrep[spin] = irr;
          }
        }
      }
    } /* end loop over irreps */

    if (!found_guess) {
      fprintf(outfile, 
              "cscf: (orb_mix): Can't find a valid guess for LUMO irrep!\n");
      exit(1);
    } 


    if (print > 1) {
      fprintf(outfile, "Identified alpha HOMO as irrep %s, energy %12.6lf\n",
              scf_info[homo_irrep[0]].irrep_label, homo_energy[0]);
      fprintf(outfile, "Identified alpha LUMO as irrep %s, energy %12.6lf\n",
              scf_info[lumo_irrep[0]].irrep_label, lumo_energy[0]);
    }


    /* now see if HOMO and LUMO have the same irrep! */
    if (homo_irrep[spin] != lumo_irrep[spin]) {
      fprintf(outfile, "Identified HOMO as irrep %s, energy %12.6lf\n",
            scf_info[homo_irrep[spin]].irrep_label, homo_energy[spin]);
      fprintf(outfile, "Identified LUMO as irrep %s, energy %12.6lf\n",
            scf_info[lumo_irrep[spin]].irrep_label, lumo_energy[spin]);
      fprintf(outfile, "\tHOMO and LUMO have different irreps.\n");
      fprintf(outfile, 
        "\tWill use lowest unoccupied from irrep %s instead\n", 
        scf_info[homo_irrep[spin]].irrep_label);
        lumo_irrep[spin] = homo_irrep[spin];
    }

    /* now do the same for beta as for alpha */
    homo_irrep[1] = homo_irrep[0];
    lumo_irrep[1] = lumo_irrep[0];

    if (print > 4) {
      for (spin=0; spin<2; spin++) {
        num_mo = scf_info[lumo_irrep[spin]].num_mo;
        fprintf(outfile,"\n C matrix before mixing (%s spin)\n",
                spin == 0 ? "alpha" : "beta");
        print_mat(spin_info[spin].scf_spin[lumo_irrep[spin]].cmat,
                  num_mo,num_mo,outfile);
      }
    }

    /* figure out mixing coefficients */
    mixing_fraction = 0.5;
    errcod = ip_data("MIXING_FRACTION","%lf",&mixing_fraction,0);
    c1 = cos(mixing_fraction*PI/2.0);
    c2 = sin(mixing_fraction*PI/2.0);
    fprintf(outfile, "        Mixing HOMO and LUMO by %d\%\n", 
            (int) (mixing_fraction * 100.0));
    fprintf(outfile, "        Mixing coefficients are %lf and %lf\n", c1, c2);
    
    /* Now mix HOMO and LUMO for each spin */
    for (spin=0; spin<2; spin++) {
      irr = homo_irrep[spin];
      num_mo = scf_info[irr].num_mo;
      homo_orb = spin_info[spin].scf_spin[irr].noccup-1;
      lumo_orb = homo_orb+1;
      for (i=0; i<num_mo; i++) {
        tval1 =  c1*spin_info[spin].scf_spin[irr].cmat[i][homo_orb]
              +  c2*spin_info[spin].scf_spin[irr].cmat[i][lumo_orb];
        tval2 = -c2*spin_info[spin].scf_spin[irr].cmat[i][homo_orb]
              +  c1*spin_info[spin].scf_spin[irr].cmat[i][lumo_orb];
        spin_info[spin].scf_spin[irr].cmat[i][homo_orb] = 
          spin == 0 ? tval1 : tval2;
        spin_info[spin].scf_spin[irr].cmat[i][lumo_orb] = 
          spin == 0 ? tval2 : tval1;
      }

      if (print > 4) {
        fprintf(outfile, "C matrix after mixing (%s spin)\n", 
                spin==0 ? "alpha" : "beta");
	print_mat(spin_info[spin].scf_spin[irr].cmat,num_mo,num_mo,outfile);
      }

    } /* end loop over spins */
    
}
   

