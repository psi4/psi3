#!/usr/bin/perl  

$TOL = 10**-5;
$REF_FILE = "output.ref";
$TEST_FILE = "output.dat";
$RESULT = "cc12.test";
$NUM_EVALS = 2;
$NUM_SYMMS = 4;

#system ("input");
#system ("psi3");

extract_data($REF_FILE,$Enuc_ref,$Ehf_ref,$Eccsd_ref);
extract_data($TEST_FILE,$Enuc_test,$Ehf_test,$Eccsd_test);

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

  seek(OUT,0,0);
  while (<OUT>) {
    if (/Total CCSD energy/) {
      @data3 = split(/ +/, $_);
      $_[3] = $data3[4];
    }
  }

  $j=0;
  $linenum=0;
  $symm = -1;
  foreach $line (@output) {
    if ($line =~ m/Symmetry of excited state:/) {
      $symm++; $eval=0;
    }
    $linenum++;
    if ($line =~ m/Largest components of/) {
      @test = split (/ +/,$output[$linenum-3]);
      if($_[0] eq $REF_FILE) {
        $evals_ref[$symm][$eval] = $test[5];
      }
      elsif($_[0] eq $TEST_FILE) {
        $evals_test[$symm][$eval] = $test[5];
      }
     $eval++;
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

  $diff_ccsd = abs ($Eccsd_ref - $Eccsd_test);
  if ($diff_ccsd > $TOL) {
    printf "\nCCSD Energy              ... FAILED\n";
  }
  else {
    printf "\nCCSD Energy              ... PASSED\n";
  }

  for($j=0; $j < $NUM_SYMMS; $j++) {
    $OK = "PASS";
    for($k=0; $k < $NUM_EVALS; $k++) {
#       printf "%15.8f %15.8f\n", $evals_ref[$j][$k], $evals_test[$j][$k];
        if(abs($evals_ref[$j][$k] - $evals_test[$j][$k]) > $TOL) {
          $OK = "FAIL";
        }
    }
    if ($OK eq "FAIL") {
      printf "\nEOM-CCSD for irrep %d ... FAILED\n", $j;
    }
    else {
      printf "\nEOM-CCSD for irrep %d ... PASSED\n", $j;
    }
  }

  close (RE);

}
