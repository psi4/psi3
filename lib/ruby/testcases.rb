# Make sure color.rb is included
require 'color'

PSITEST_ETOL = 10e-8;           # Default test criterion for energies
PSITEST_ENUCTOL = 10e-10;       # Check nuclear repulsion energy tighter than other energies

def test_scf_energy(expected)
  retcode = 0
  
  # Make sure SCF is in the checkpoint
  if Psi::Chkpt::exist?("SCF energy")
    if (Psi::Chkpt::escf - expected).abs < PSITEST_ETOL
      puts "SCF final energy:          " + green("PASSED")
    else
      puts "SCF final energy:          " + red("FAILED")
      puts "Obtained: #{Psi::Chkpt::escf}"
      puts "Expected: #{expected}"
      retcode = 1
    end
  else
    puts red("Error: ") + "SCF energy not found in checkpoint."
    retcode = 1
  end
  
  return retcode
end

def test_ccsd_energy(expected)
  retcode = 0

  # Make sure CCSD is in the checkpoint
  if Psi::Chkpt::exist?("CCSD Energy")
    if (Psi::Chkpt::eccsd - expected).abs < PSITEST_ETOL
      puts "CCSD final energy:         " + green("PASSED")
    else
      puts "CCSD final energy:         " + red("FAILED")
      puts "Obtained: #{Psi::Chkpt::eccsd}"
      puts "Expected: #{expected}"
      retcode = 1
    end
  else
    puts red("Error: ") + "CCSD energy not found in checkpoint."
    retcode = 1
  end

  return retcode
end

def test_t_energy(expected)
  retcode = 0

  # Make sure (T) is in the checkpoint
  if Psi::Chkpt::exist?("(T) Energy")
    if (Psi::Chkpt::e_t - expected).abs < PSITEST_ETOL
      puts "(T) final energy:          " + green("PASSED")
    else
      puts "(T) final energy:          " + red("FAILED")
      puts "Obtained: #{Psi::Chkpt::e_t}"
      puts "Expected: #{expected}"
      retcode = 1
    end
  else
    puts red("Error: ") + "(T) energy not found in checkpoint."
    retcode = 1
  end

  return retcode
end

def test_nuclear_repulsion(expected)
  retcode = 0
  
  # Make sure it is found in the checkpoint
  if Psi::Chkpt::exist?("Nuclear rep. energy")
    if (Psi::Chkpt::enuc - expected).abs < PSITEST_ENUCTOL
      puts "Nuclear repulsion energy:  " + green("PASSED")
    else
      puts "Nuclear repulsion energy:  " + red("FAILED")
      puts "Obtained: #{Psi::Chkpt::enuc}"
      puts "Expected: #{expected}"
      retcode = 1
    end
  else
    puts red("Error:") + "Nuclear repulsion energy not found in checkpoint."
    retcode = 1
  end
  
  return retcode
end

def test_total_energy(expected)
  retcode = 0
  
  # Make sure it is found in the checkpoint
  if Psi::Chkpt::exist?("Total energy")
    if (Psi::Chkpt::etot - expected).abs < PSITEST_ETOL
      puts "Total energy:              " + green("PASSED")
    else
      puts "Total energy:              " + red("FAILED")
      puts "Obtained: #{Psi::Chkpt::etot}"
      puts "Expected: #{expected}"
      retcode = 1
    end
  else
    puts red("Error:") + "Total energy not found in checkpoint."
    retcode = 1
  end
  
  return retcode
end
