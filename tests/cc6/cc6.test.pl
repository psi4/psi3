#!/usr/bin/perl  

$TOL = 10**-8;
$REF_FILE = "output.ref";
$TEST_FILE = "output.dat";
$RESULT = "cc6.test";

system ("input");
system ("psi3");

extract_data($REF_FILE,$Enuc_ref,$Ehf_ref,$Eccsd_ref,$Eccsd_t_ref);
extract_data($TEST_FILE,$Enuc_test,$Ehf_test,$Eccsd_test,$Eccsd_t_test);

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

  seek(OUT,0,0);
  while (<OUT>) {
    if (/Total CCSD[\(]T[\)] energy/) {
      @data4 = split(/ +/, $_);
      $_[4] = $data4[4];
    }
  }

  close(OUT);
}
###############################################################################
sub compare_data
{
  open(OUT, ">$RESULT") || die "cannot open $RESULT: $!";

  select (RE);

  printf "\nCC6:\n";

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

  $diff_ccsd_t = abs ($Eccsd_t_ref - $Eccsd_t_test);
  if ($diff_ccsd_t > $TOL) {
    printf "\nCCSD(T) Energy           ... FAILED\n\n";
  }
  else {
    printf "\nCCSD(T) Energy           ... PASSED\n\n";
  }

  close (RE);
}
