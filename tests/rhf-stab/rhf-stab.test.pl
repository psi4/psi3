#!/usr/bin/perl  

$TOL = 10**-8;
$STAB_TOL = 10**-4;
$REF_FILE = "output.ref";
$TEST_FILE = "output.dat";
$RESULT = "rhf-stab.test";
$NUM_EVALS = 5;
$NUM_SYMMS = 4;

system ("input");
system ("psi3");

extract_data($REF_FILE,$Enuc_ref,$Escf_ref);
extract_data($TEST_FILE,$Enuc_test,$Escf_test);

compare_data();

###############################################################################
sub extract_data
{
  open(OUT, "$_[0]") || die "cannot open $_[0] $!";

  seek(OUT,0,0);
  while(<OUT>) {
    if(/Nuclear Repulsion Energy    =/) {
      @data1 = split(/ +/, $_);
      $_[1] = $data1[4];
    }
  }

  seek(OUT,0,0);
  while (<OUT>) {
    if (/SCF total energy       =/) {
      @data2 = split(/ +/, $_);
      $_[2] = $data2[4];
    }
  }
  close (OUT);


# scan the output (stored as an array) for the stability eigenvalues
  $linenum = 0;
  $start1 = 0;
  $start2 = 0;

  open(OUT, "$_[0]") || die "cannot open $_[0] $!";
  @infile = <OUT>;
  close (OUT);

  for $line (@infile) {
    $linenum++;
    if ($line =~ m/RHF->RHF/) {
      $start1 = $linenum;
    }
    if ($line =~ m/RHF->UHF/) {
      $start2 = $linenum;
    }
  }

  for($j=0; $j<$NUM_EVALS; $j++) {
    @line = split (/ +/, $infile[$start1 + 3 + $j]);
    if($_[0] eq $REF_FILE) {
      for($k=0; $k<$NUM_SYMMS; $k++) {
        $stab1_ref[$j][$k] = $line[$k+2];
      }
    }
    elsif($_[0] eq $TEST_FILE) {
      for($k=0; $k<$NUM_SYMMS; $k++) {
        $stab1_test[$j][$k] = $line[$k+2];
      }
    }
  }

  for($j=0; $j<$NUM_EVALS; $j++) {
    @line = split (/ +/, $infile[$start2 + 3 + $j]);
    if($_[0] eq $REF_FILE) {
      for($k=0; $k<$NUM_SYMMS;$k++) {
        $stab2_ref[$j][$k] = $line[$k+2];
      }
    }
    elsif($_[0] eq $TEST_FILE) {
      for($k=0; $k<$NUM_SYMMS;$k++) {
        $stab2_test[$j][$k] = $line[$k+2];
      }
    }
  }


}
###############################################################################
sub compare_data
{
  open(RE, ">$RESULT") || die "cannot open $RESULT $!";
  select(RE);
   
  printf "\nRHF-STAB:\n";

  $diff_enuc = abs ($Enuc_ref - $Enuc_test);
  if ($diff_enuc > $TOL) {
    printf "\nNuclear Repulsion Energy  ... FAILED\n";
  }
  else {
    printf "\nNuclear Repulsion Energy  ... PASSED\n";
  }

  $diff_scf = abs ($Escf_ref - $Escf_test);
  if ($diff_scf > $TOL) {
    printf "\nSCF Energy           ... FAILED\n";
  }
  else {
    printf "\nSCF Energy           ... PASSED\n";
  }

  $stab1_out = "PASS";
  $stab2_out = "PASS";
  for($j=0; $j < $NUM_EVALS; $j++) {
    for($k=0; $k < $NUM_SYMMS; $k++) {
#       printf "%6.4f  %6.4f \n",$stab1_ref[$j][$k], $stab1_test[$j][$k];
       if(abs($stab1_ref[$j][$k] - $stab1_test[$j][$k]) > $STAB_TOL) {
          $stab1_out = "FAIL";
       }
#       printf "%6.4f  %6.4f \n",$stab2_ref[$j][$k], $stab2_test[$j][$k];
       if(abs($stab2_ref[$j][$k] - $stab2_test[$j][$k]) > $STAB_TOL) {
          $stab2_out = "FAIL";
       }
    }
  }

  if ($stab1_out eq "FAIL") {
    printf "\nRHF->RHF stability   ... FAILED\n";
  }
  else {
    printf "\nRHF->RHF stability   ... PASSED\n";
  }

  if ($stab2_out eq "FAIL") {
    printf "\nRHF->UHF stability   ... FAILED\n";
  }
  else {
    printf "\nRHF->UHF stability   ... PASSED\n";
  }

  close (RE);
}
###############################################################################
