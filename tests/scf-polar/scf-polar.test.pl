#!/usr/bin/perl  

$TOL = 10**-10;
$PTOL = 10**-4;
$REF_FILE = "output.ref";
$TEST_FILE = "output.dat";
$NDOF = 3;
$RESULT = "scf-polar.test";

system ("input");
system ("psi3");

extract_data($REF_FILE,$Enuc_ref,$Ehf_ref,$Eref_ref);
extract_data($TEST_FILE,$Enuc_test,$Ehf_test,$Eref_test);

compare_data();

###############################################################################
sub extract_data
{
  open(OUT, "$_[0]") || die "cannot open $_[0]: $!";
  
  @infile = <OUT>;

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
    if (/Reference energy             =/) {
      @data3 = split(/ +/, $_);
      $_[3] = $data3[3];
    }
  }
  
  close(OUT);

  $j=0;
  $linenum=0;
  foreach $line (@infile) {
    $linenum++;
    if ($line =~ m/Polarizability Tensor/) {
      while ($j<3) {
        @test = split (/ +/,$infile[$linenum+4+$j]);
        if ($_[0] eq $REF_FILE) {
          $a_ref[$j] = $test[2];
          $b_ref[$j] = $test[3];
          $c_ref[$j] = $test[4];
        }
        elsif ($_[0] eq $TEST_FILE) {
          $a_test[$j] = $test[2];
          $b_test[$j] = $test[3];
          $c_test[$j] = $test[4];
        }
        $j++;
      }
    }
  }
 
}
###############################################################################
sub compare_data
{
  open (RE, ">$RESULT") || die "cannot open $RESULT: $!";

  select (RE);

  printf("Enuc(ref)  = %20.10f\n", $Enuc_ref);
  printf("Enuc(test) = %20.10f\n", $Enuc_test);
  printf("Escf(ref)  = %20.10f\n", $Ehf_ref);
  printf("Escf(test) = %20.10f\n", $Ehf_test);
  printf("Eref(ref)  = %20.10f\n", $Eref_ref);
  printf("Eref(test) = %20.10f\n", $Eref_test);

  printf "\nSCF-POLAR:\n";
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

  $diff_ref = abs ($Eref_ref - $Eref_test);
  if ($diff_ref > $TOL) {
    printf "\nReference Energy         ... FAILED\n";
  }
  else {
    printf "\nReference Energy         ... PASSED\n";
  }


  for ($i=0; $i<$NDOF; $i++) {
    $diff_a[$i] = abs ($a_ref[$i] - $a_test[$i]);
    $diff_b[$i] = abs ($b_ref[$i] - $b_test[$i]);
    $diff_c[$i] = abs ($c_ref[$i] - $c_test[$i]);
  }

  for ($i=0; $i<$NDOF; $i++) {
    if ($diff_a[$i] > $PTOL) {
      $pout = "FAIL";
    }
    if ($diff_b[$i] > $PTOL) {
      $pout = "FAIL";
    }
    if ($diff_c[$i] > $PTOL) {
      $pout = "FAIL";
    }
  }
  
  if ($pout eq "FAIL") {
    printf "\nPolarizability           ... FAILED\n";
  }
  else {
    printf "\nPolarizability           ... PASSED\n";
  } 
}
###############################################################################
