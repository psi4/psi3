#!/usr/bin/perl  

require("../psitest.pl");

$TOL = 10**-8;
$EVAL_TOL = 10**-9;
$REF_FILE = "output.ref";
$TEST_FILE = "output.dat";
$RESULT = "cc17.test";
$NUM_EVALS = 2;
$NUM_SYMMS = 4;

system ("input");
system ("psi3");

open(RE, ">$RESULT") || die "cannot open $RESULT: $!";
select (RE);
printf "\nCC17:\n";

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

if (abs(seek_ccsd($REF_FILE) - seek_ccsd($TEST_FILE)) > $TOL) {
  fail_test("CCSD Energy");
}
else {
  pass_test("CCSD Energy");
}

@eom_ref = seek_eomcc($REF_FILE);
@eom_test = seek_eomcc($TEST_FILE);

if(!compare_arrays(\@eom_ref,\@eom_test,$NUM_SYMMS,$NUM_EVALS,$EVAL_TOL)) {
  fail_test("EOM-CCSD energy");
}
else {
  pass_test("EOM-CCSD energy");
}
close (RE);
