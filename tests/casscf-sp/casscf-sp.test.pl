#!/usr/bin/perl  

$TOL = 10**-10;
$REF_FILE = "file14.ref";
$TEST_FILE = "file14.dat";
$REF_OUT = "output.ref";
$TEST_OUT = "output.dat";
$RESULT = "casscf-sp.test";

system ("input");
system ("psi3");

extract_data($REF_OUT,$REF_FILE,$Enuc_ref,$Ehf_ref,$Ecas_ref);
extract_data($TEST_OUT,$TEST_FILE,$Enuc_test,$Ehf_test,$Ecas_test);

compare_data();

###############################################################################
sub extract_data
{
  open(OUT, "$_[0]") || die "cannot open $_[0]: $!";

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
  
  close (OUT);

  open(OUT, "$_[1]") || die "cannot open $_[1]: $!";
  
  @cas = <OUT>;

  @line1 = split(/ +/, $cas[0]);
  $niter = $line1[1];
  @line2 = split(/ +/, $cas[$niter]);
  $_[4] = $line2[5];
  close(OUT);

}
###############################################################################
sub compare_data
{
  open(RE, ">$RESULT") || die "cannot open $RESULT: $!";

  select (RE);

  printf "\nCASSCF-SP:\n";
  
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

  $diff_cas = abs ($Ecas_ref - $Ecas_test);
  if ($diff_cas > $TOL) {
    printf "\nCASSCF Energy            ... FAILED\n\n";
  }
  else {
    printf "\nCASSCF Energy            ... PASSED\n\n";
  }

  close (RE);
}
###############################################################################
