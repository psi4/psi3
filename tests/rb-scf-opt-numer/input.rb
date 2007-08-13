# SCF geometry optimization by energy points
#  This test case includes extra commands to do a complete test

# The values to be expected
ETOTAL = -76.0236150188265
ENUCTOTAL = 9.32931478014697

# Make sure we are running in PSI
if $0 != "Psi"
  puts "Error: Not running within Psi."
  puts "       This file must be run wih: psirb inputfilename"
  exit 1
end

# Attempt to load psi3 and testcases
begin
  require 'psi3'
rescue LoadError => e
  puts "Error: Unable to load psi3.rb."
  puts "       Make sure that either the binaries are installed"
  puts "       in the correct place, or that PSIDATADIR is set."
  exit 1
end

# Attempt to load psi3 and testcases
begin
  require 'testcases'
rescue LoadError => e
  puts "Error: Unable to load testcases.rb."
  puts "       Make sure that either the binaries are installed"
  puts "       in the correct place, or that PSIDATADIR is set."
  exit 1
end

# Tell Psi not to print
Psi::quiet = true

# Everything we need is loaded.
# Do the equivalient of the scf-opt-numer test
Psi::label = "6-31G** SCF H2O"
Psi::wavefunction = "scf"
Psi::reference = "rhf"
Psi::basis = "6-31G**"
Psi::zmat = [
  [ O ],
  [ H, 1, 1.00 ],
  [ H, 1, 1.00, 2, 103.1 ]
]

# Do the optimization
optking(:gradients => false) do
  rhf
end

# Clean out the scratch space
clean

Psi::quiet = false

retcode = 0
retcode |= test_scf_energy(ETOTAL)
retcode |= test_nuclear_repulsion(ENUCTOTAL)

exit retcode
