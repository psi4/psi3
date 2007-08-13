# Reproduces test case: fci-h2o
#  This test case includes extra commands to do a complete test
ESCF = -75.985323665263
ETOTAL = -76.1210978474779
ENUC = 9.234218520912

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
  [ O,       0.0000000000,        0.0000000000,       -0.0742719254],
  [ H,       0.0000000000,       -1.4949589982,       -1.0728640373],
  [ H,       0.0000000000,        1.4949589982,       -1.0728640373]
]

# Set the basis
Psi::basis = "6-31G"
Psi::reference = "rhf"
Psi::wavefunction = "detci"

# Clean Psi scratch files
clean

# Run input
input("units" => "bohr")

# Run scf (which internally runs cints)
scf("docc" => [3, 0, 1, 1])

detci("fci" => "true")

# Clean the scratch out
clean

Psi::quiet = false

retcode = 0
retcode |= test_nuclear_repulsion(ENUC)
retcode |= test_scf_energy(ESCF)
retcode |= test_total_energy(ETOTAL)

exit retcode
