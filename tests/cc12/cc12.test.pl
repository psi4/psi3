#!/usr/bin/perl  

$TOL = 10**-8;
$REF_FILE = "output.ref";
$TEST_FILE = "output.dat";
$RESULT = "cc12.test";

system ("input");
system ("psi3");

extract_data($REF_FILE,$Enuc_ref,$Ehf_ref,$Eeom_a1_ref,
             $Eeom_a2_ref,$Eeom_b1_ref,$Eeom_b2_ref);
extract_data($TEST_FILE,$Enuc_test,$Ehf_test,$Eeom_a1_test,
             $Eeom_a2_test,$Eeom_b1_test,$Eeom_b2_test);

compare_data();

###############################################################################
sub extract_data
{
  open(OUT, "$_[0]") || die "cannot open $_[0]: $!";

  @output = <OUT>;

  seek(OUT,0,0);
  while (<OUT>) {
    if (/Nuclear Repulsion Energy    =/) {
      @data1 = split(/ +/, $_);
      $_[1] = $data1[4];
    }
  }

  seek(OUT,0,0);
  while (<OUT>) {
    if (/SCF total energy   =/) {
      @data2 = split(/ +/, $_);
      $_[2] = $data2[4];
    }
  }

  $j=0;
  $linenum=0;
  foreach $line (@output) {
    $linenum++;
    if ($line =~ m/Final Energetic Summary/) {
      @test = split (/ +/,$output[$linenum+2]);
      $_[$j+3] = $test[5];
      $j++;
    }
  }

  close (OUT);
}
###############################################################################
sub compare_data
{
  open(RE, ">$RESULT") || die "cannot open $RESULT: $!";

  select (RE);

  printf "\nCC12:\n";

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

  $diff_eom_a1 = abs ($Eeom_a1_ref - $Eeom_a1_test);
  if ($diff_eom_a1 > $TOL) {
    printf "\nEOM-CCSD A1 Energy       ... FAILED\n";
  }
  else {
    printf "\nEOM-CCSD A1 Energy       ... PASSED\n";
  }

  $diff_eom_a2 = abs ($Eeom_a2_ref - $Eeom_a2_test);
  if ($diff_eom_a2 > $TOL) {
    printf "\nEOM-CCSD A2 Energy       ... FAILED\n";
  }
  else {
    printf "\nEOM-CCSD A2 Energy       ... PASSED\n";
  }

  $diff_eom_b1 = abs ($Eeom_b1_ref - $Eeom_b1_test);
  if ($diff_eom_b1 > $TOL) {
    printf "\nEOM-CCSD B1 Energy       ... FAILED\n";
  }
  else {
    printf "\nEOM-CCSD B1 Energy       ... PASSED\n";
  }

  $diff_eom_b2 = abs ($Eeom_b2_ref - $Eeom_b2_test);
  if ($diff_eom_b2 > $TOL) {
    printf "\nEOM-CCSD B2 Energy       ... FAILED\n";
  }
  else {
    printf "\nEOM-CCSD B2 Energy       ... PASSED\n";
  }

  close (RE);

}
