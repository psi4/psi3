#!/usr/bin/perl  

$TOL = 10**-10;
$REF_FILE = "file11.ref";
$TEST_FILE = "psi.file11.dat";
$RESULT = "scf-opt.test";

system ("input");
system ("psi3");

extract_data($REF_FILE,$Escf_ref);
extract_data($TEST_FILE,$Escf_test);

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
  $_[1] = $line[2];

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
}
###############################################################################
sub compare_data
{
  open(RE, ">$RESULT") || die "cannot open $RESULT $!";
  select(RE);
   
  printf "\nSCF-OPT:\n";

  $diff_scf = abs ($Escf_ref - $Escf_test);
  if ($diff_ccsd > $TOL) {
    printf "\nSCF Energy           ... FAILED\n";
  }
  else {
    printf "\nSCF Energy           ... PASSED\n";
  }

  for ($i=0; $i<$natom; $i++) {
    $diff_x[$i] = abs ($x_ref[$i] - $x_test[$i]);
    $diff_y[$i] = abs ($y_ref[$i] - $y_test[$i]);
    $diff_z[$i] = abs ($z_ref[$i] - $z_test[$i]);
  }

  for ($i=0; $i<$natom; $i++) {
    if ($diff_x[$i] > $TOL) {
      $out = "FAIL";
    }
    if ($diff_y[$i] > $TOL) {
      $out = "FAIL";
    }
    if ($diff_z[$i] > $TOL) {
      $out = "FAIL";
    }
  }

  if ($out eq "FAIL") {
    printf "\nGeometry             ... FAILED\n\n";
  }
  else {
    printf "\nGeometry             ... PASSED\n\n";
  }

  close (RE);
}
###############################################################################
