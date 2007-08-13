# CCSD single point
#  This test case includes extra commands to do a complete test.

ESCF = -1.100153777808
ECCSD = -0.039919706944083
ENUC = 0.5291772490000

# Make sure we are running in Psi
if $0 != "Psi"
  puts "Not running within Psi!!!"
  exit
end

# Attempt to load psi3
begin
  require 'psi3'
rescue LoadError => e
  puts "Error: Unable to load psi3.rb"
  puts "       Make sure that either the binaries are installed"
  puts "       in the correct place, or that PSIDATADIR is set."
  exit 1
end

# Attempt to load testcases
begin
  require 'testcases'
rescue LoadError => e
  puts "Error: Unable to load testcases.rb"
  puts "       Make sure that either the binaries are installed"
  puts "       in the correct place, or that PSIDATADIR is set."
  exit 1
end

Psi::quiet = true

# Set the geometry
Psi::geometry = [
  [ "H", 0.0, 0.0, 0.0 ],
  [ "H", 0.0, 0.0, 1.0 ]
]

# Set the basis
Psi::basis = "cc-pVDZ"

# Clean Psi scratch files
clean

# Set the final desired wavefunction
Psi::wavefunction = "ccsd_t"

# Run input
input

# Run scf (which internally runs cints)
rhf

# Run ccenergy to compute ccsd energy
ccsd

# Clean the scratch out
clean

Psi::quiet = false

retcode = 0
retcode |= test_nuclear_repulsion(ENUC)
retcode |= test_scf_energy(ESCF)
retcode |= test_ccsd_energy(ECCSD)

exit retcode

