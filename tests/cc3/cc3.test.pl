#!/usr/bin/perl  

$TOL = 10**-8;
$GTOL = 10**-8;
$HTOL = 10**-1;
$REF_FILE = "file11.ref";
$TEST_FILE = "psi.file11.dat";
$REF_OUT = "output.ref";
$TEST_OUT = "output.dat";
$RESULT = "cc3.test";
$NDOF = 3;

system ("input");
system ("psi3");

extract_data($REF_FILE,$REF_OUT,$Eccsd_ref);
extract_data($TEST_FILE,$TEST_OUT,$Eccsd_test);

compare_data();

###############################################################################
sub extract_data
{
  # Open input file, place into array for parsing, close input file
  open(IN, "$_[0]") || die "cannot open $_[0] $!";
  @infile = <IN>; #file11 is copied to @outfile
  close (IN);

  # Initialize variables
  $match = "SCF";
  $linenum = 0;
  $lastiter = 0;
  $flag = 0;

  foreach $line (@infile) {
    $linenum++;
    if ($line =~ m/$match/) {
      if ($flag == 0) {
        $flag = 1;
      }
      $lastiter = $linenum;
    }
  }

  @line = split (/ +/, $infile[$lastiter]);
  $natom = $line[1];
  $_[2] = $line[2];

  $j=1;
  while ($j<$natom+1) {
    @line2 = split (/ +/, $infile[$lastiter+$j]);
    if ($_[0] eq $REF_FILE) {
      $x_ref[$j-1] = $line2[2];
      $y_ref[$j-1] = $line2[3];
      $z_ref[$j-1] = $line2[4];
    }
    elsif ($_[0] eq $TEST_FILE) {
      $x_test[$j-1] = $line2[2];
      $y_test[$j-1] = $line2[3];
      $z_test[$j-1] = $line2[4];
    }
    $j++;
  }

  open(OUT, "$_[1]") || die "cannot open $_[1]: $!";

  @infile = <OUT>;

  close(OUT);

  $j=0;
  $linenum=0;
  foreach $line (@infile) {
    $linenum++;
    if ($line =~ m/Harmonic Vibrational Frequencies/) {
      while ($j<$NDOF) {
        @test = split (/ +/,$infile[$linenum+$j]);
        if ($_[0] eq $REF_FILE) {
          $hvf_ref[$j] = $test[2];
        }
        elsif ($_[0] eq $TEST_FILE) {
          $hvf_test[$j] = $test[2];
        }
        $j++;
      }
    }
  }

}
###############################################################################
sub compare_data
{
  open(RE, ">$RESULT") || die "cannot open $RESULT: $!";

  select (RE);

  printf "\nCC3:\n";

  $diff_ccsd = abs ($Eccsd_ref - $Eccsd_test);

  if ($diff_ccsd > $TOL) {
    printf "\nCCSD Energy                            ... FAILED\n";
  }
  else {
    printf "\nCCSD Energy                            ... PASSED\n";
  }

  for ($i=0; $i<$natom; $i++) {
    $diff_x[$i] = abs ($x_ref[$i] - $x_test[$i]);
    $diff_y[$i] = abs ($y_ref[$i] - $y_test[$i]);
    $diff_z[$i] = abs ($z_ref[$i] - $z_test[$i]);
  }

  for ($i=0; $i<$natom; $i++) {
    if ($diff_x[$i] > $GTOL) {
      $out = "FAIL";
    }
    if ($diff_y[$i] > $GTOL) {
      $out = "FAIL";
    }
    if ($diff_z[$i] > $GTOL) {
      $out = "FAIL";
    }
  }

  if ($out eq "FAIL") {
    printf "\nGeometry                               ... FAILED\n";
  }
  else {
    printf "\nGeometry                               ... PASSED\n";
  }

  for ($i=0; $i<$NDOF; $i++) {
    $diff_hvf[$i] = abs ($hvf_ref[$i] - $hvf_test[$i]);
  }

  for ($i=0; $i<$NDOF; $i++) {
    if ($diff_hvf[$i] > $HTOL) {
      $hout = "FAIL";
    }
  }

  if ($hout eq "FAIL") {
    printf "\nHarmonic Vibrational Frequencies       ... FAILED\n\n";
  }
  else {
    printf "\nHarmonic Vibrational Frequencies       ... PASSED\n\n";
  }

  close (RE);
}
###############################################################################
