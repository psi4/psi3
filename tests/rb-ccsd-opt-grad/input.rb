# CCSD geometry optimization by analytic gradients
#  This test case includes extra commands to do a complete test

# The values to be expected
ESCF = -76.0229402438764
ECCSD = -0.208238239222441
ETOTAL = -76.231178483098901
ENUCTOTAL = 9.16617496476892

# Make sure we are running in PSI
if $0 != "Psi"
  puts "Error: Not running within Psi."
  puts "       This file must be run wih: psirb inputfilename"
  exit 1
end

# Attempt to load psi3
begin
  require 'psi3'
rescue LoadError => e
  puts "Error: Unable to load psi3.rb."
  puts "       Make sure that either the binaries are installed"
  puts "       in the correct place, or that PSIDATADIR is set."
  exit 1
end

# Attempt to load testcases
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
Psi::label = "6-31G** CCSD H2O"
Psi::wavefunction = "ccsd"
Psi::reference = "rhf"
Psi::basis = "6-31G**"
Psi::zmat = [
  [ O ],
  [ H, 1, 1.00 ],
  [ H, 1, 1.00, 2, 103.1 ]
]

# Do the optimization
optking do
  rhf
  ccsd
end

# Clean out the scratch space
clean

Psi::quiet = false

retcode = 0
retcode |= test_nuclear_repulsion(ENUCTOTAL)
retcode |= test_scf_energy(ESCF)
retcode |= test_ccsd_energy(ECCSD)
retcode |= test_total_energy(ETOTAL)

exit retcode
