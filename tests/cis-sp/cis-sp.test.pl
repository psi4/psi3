#!/usr/bin/perl  

require("../psitest.pl");

$TOL = 10**-8;
$EVAL_TOL = 10**-9;
$REF_FILE = "output.ref";
$TEST_FILE = "output.dat";
$RESULT = "cis-sp.test";
$NUM_EVALS = 2;
$NUM_SYMMS = 4;

system ("psi3");

open(RE, ">$RESULT") || die "cannot open $RESULT: $!";
select (RE);
printf "\nCIS-SP:\n";

if (abs(seek_nuc($REF_FILE) - seek_nuc($TEST_FILE)) > $TOL) {
  fail_test("Nuclear Repulsion Energy");
}
else {
  pass_test("Nuclear Repulsion Energy");
}

if (abs(seek_scf($REF_FILE) - seek_scf($TEST_FILE)) > $TOL) {
  fail_test("SCF Energy");
}
else {
  pass_test("SCF Energy");
}

@cis_ref = seek_cis($REF_FILE);
@cis_test = seek_cis($TEST_FILE);

if(!compare_arrays(\@cis_ref,\@cis_test,$NUM_SYMMS,$NUM_EVALS,$EVAL_TOL)) {
  fail_test("CIS Energies");
}
else {
  pass_test("CIS Energies");
}
close (RE);
