** RAK, Feb. 2005, notes on parallel coarse graining
** NJD, Feb. 2005, added some practical considerations

Here is a recipe for doing gradients by energy points, allowing the
energies to be computed separately and coarse-grain parallelized.
In order to make sure your checkpoint and optking files do not get
lost or overwritten, it is best to save these files in your work directory
by making sure the following is in your "input.dat" file:

  files : (
    file32: ( nvolume = 1  volume1 = "./" )
    file1: ( nvolume = 1  volume1 = "./" )
  )

This example is for SCF. See notes below for other types of
finite differences.

**** Generate the displacements in the root directory,
input
cints
cscf

**** and run one of the following

optking --disp_irrep --irrep 1    ** for optimizations,
optking --disp_freq_energy_cart   ** for frequencies from energies, or
optking --disp_freq_grad_cart     ** for frequencies from gradients

**** Repeat the following for all i displacements:
  optking --disp_load
  input --keepchkpt --chkptgeom --noreorient
  
  **make a directory for displacement i
      Within the displacement directory: {
          **copy the root checkpoint file into the displacement directory
          **copy the root input file into the displacement directory
          **modify the input file so that only an energy is calculated, beginning
            with the integral calculation, e.g. exec = ("cints" "cscf")
      }

  ** In the root directory, increment the current displacement number by running
  optking --disp_num_plus 
  (You don't need to do this after the final displacement has been made.)
****


**** All the energy calculations can then be run parallel ****
When done, put all energies into a text file named "energy.dat" - one energy on each line.
Do not include the reference energy in "energy.dat"

**** Compute the gradient ****
** Reset the checkpoint prefix by running
  optking --reset_prefix
** Compute the gradient with **
  optking --grad_energy --energy_dat
** Take optimization step with **
  optking --opt_step


*****
To get frequencies from gradients, follow a similar procedure,
but the last step is 

optking --freq_grad_cart --grad_dat

which tells optking to read the gradients from a file11.dat
file that the user will build from the displaced gradients.



*****
To get frequencies from energies, follow a similar procedure,
but the last step is

optking --freq_energy_cart --energy_dat
