#csh
# 
# This file contains a script that illustrates coarse-grain parallelization of
# frequencies by energy points.  It requires the following:
#   1) an input.dat which will compute the undisplaced energy and the
#   displaced geometries.  For SCF, this is
#     exec = ( "input" "cints" "cscf" "optking --disp_freq_energy_cart")
#
#   2) an input.dat.disp which is an input file that will compute the energies
#   for the displacements.  This file should NOT call "input", which is done
#   by the script below.  For SCF, the file contains exec = ( "cints" " cscf").
#
#   3) an awk script that will grab the desired energy from output.dat
#   In my humble awk and for SCF this takes the form
#       BEGIN { scf_energy = 0.0; }
#       {
#         if($1 == "SCF" && $2 == "total" && $3 == "energy" ) {
#           printf("%s\n", $5);
#         }
#       }
# ** Good luck!  RAK, Feb. 2005

#script for testing parallel frequencies by energies
set ndisp=53

#run energy and make displacements
psi3

#make directories for displacements
set i = 1
while ($i <= $ndisp)

  optking --disp_load
  input --keepchkpt --chkptgeom --noreorient

  set dirname = disp$i
  mkdir $dirname

  cd $dirname
  cp ../psi.32 .
  cp ../input.dat.disp ./input.dat
  cd ..

  if ($i != $ndisp) then
    optking --disp_num_plus
  endif
  @ i ++
end

#run PSI3
set i = 1
while ($i <= $ndisp)
  echo "Running displacement" $i
  set dirname = disp$i
  cd $dirname
  psi3
  cd ..
  @ i ++
end

#collect energies
set i = 1
while ($i <= $ndisp)
  set dirname = disp$i
#the awk script looks for "SCF total energy"
  awk -f ~/lib/awk/scf_energy $dirname/output.dat >> energy.dat
  @ i ++
end

set i = 1
while ($i <= $ndisp)
  set dirname = disp$i
  rm -rf $dirname
  @ i ++
end

optking --reset_prefix
optking --freq_energy_cart --energy_dat

