#!/usr/bin/perl  

$TOL = 10**-10;
$REF_FILE = "output.ref";
$TEST_FILE = "output.dat";
$RESULT = "cisd-sp.test";

system ("input");
system ("psi3");

extract_data($REF_FILE,$Enuc_ref,$Ehf_ref,$Eci_ref);
extract_data($TEST_FILE,$Enuc_test,$Ehf_test,$Eci_test);

compare_data();

###############################################################################
sub extract_data
{
  open(OUT, "$_[0]") || die "cannot open $_[0]: $!";

  seek(OUT,0,0);
  while (<OUT>) {
    if (/Nuclear Repulsion Energy    =/) {
      @data1 = split(/ +/, $_);
      $_[1] = $data1[4];
    }
  }

  seek(OUT,0,0);
  while (<OUT>) {
    if (/total energy       =/) {
      @data2 = split(/ +/, $_);
      $_[2] = $data2[4];
    }
  }

  seek(OUT,0,0);
  while (<OUT>) {
    if (/ROOT 1/) {
      @data3 = split(/ +/, $_);
      $_[3] = $data3[4];
    }
  }

  close(OUT);
}
###############################################################################
sub compare_data
{
  open(RE, ">$RESULT") || die "cannot open $RESULT: $!";

  select (RE);

  printf "\nCISD-SP:\n";
  
  $diff_nuc = abs ($Enuc_ref - $Enuc_test);
  if ($diff_nuc > $TOL) {
    printf "\nNuclear Repulsion Energy ... FAILED\n";
  }
  else {
    printf "\nNuclear Repulsion Energy ... PASSED\n";
  }

  $diff_hf = abs ($Ehf_ref - $Ehf_test);
  if ($diff_hf > $TOL) {
    printf "\nRHF Energy               ... FAILED\n";
  }
  else {
    printf "\nRHF Energy               ... PASSED\n";
  }

  $diff_ci = abs ($Eci_ref - $Eci_test);
  if ($diff_ci > $TOL) {
    printf "\nCISD Energy              ... FAILED\n\n";
  }
  else {
    printf "\nCISD Energy              ... PASSED\n\n";
  }

  close (RE);
}
