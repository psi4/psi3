#!/usr/bin/perl  

$TOL = 10**-10;
$HTOL = 10**-2;
$ITOL = 10**-3;
$REF_FILE = "output.ref";
$TEST_FILE = "output.dat";
$NDOF = 9;
$RESULT = "scf-freq.test";

system ("input");
system ("psi3");

extract_data($REF_FILE,$Enuc_ref,$Ehf_ref);
extract_data($TEST_FILE,$Enuc_test,$Ehf_test);

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
      $_[2] = $data1[4];
    }
  }

  seek(OUT,0,0);
  while (<OUT>) {
    if (/total energy       =/) {
      @data2 = split(/ +/, $_);
      $_[3] = $data2[4];
    }
  }
  
  close(OUT);

  $j=0;
  $linenum=0;
  foreach $line (@infile) {
    $linenum++;
    if ($line =~ m/Harmonic Frequency/) {
      while ($j<$NDOF) {
        @test = split (/ +/,$infile[$linenum+2+$j]);
        if ($_[0] eq $REF_FILE) {
          $hvf_ref[$j] = $test[2];
          $int_ref[$j] = $test[3];
        }
        elsif ($_[0] eq $TEST_FILE) {
          $hvf_test[$j] = $test[2];
          $int_test[$j] = $test[3];
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

  printf "\nSCF-FREQ:\n";
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

  for ($i=0; $i<$NDOF; $i++) {
    $diff_hvf[$i] = abs ($hvf_ref[$i] - $hvf_test[$i]);
    $diff_int[$i] = abs ($int_ref[$i] - $int_test[$i]);
  }

  for ($i=0; $i<$NDOF; $i++) {
    if ($diff_hvf[$i] > $HTOL) {
      $hout = "FAIL";
    }
    if ($diff_int[$i] > $ITOL) {
      $iout = "FAIL";
    }
  }
  
  if ($hout eq "FAIL") {
    printf "\nHarmonic Frequency       ... FAILED\n";
  }
  else {
    printf "\nHarmonic Frequency       ... PASSED\n";
  } 

  if ($iout eq "FAIL") {
    printf "\nInfrared Intensity       ... FAILED\n\n";
  }
  else {
    printf "\nInfrared Intensity       ... PASSED\n\n";
  }

}
###############################################################################
