#!/usr/bin/perl  

require("../psitest.pl");

$TOL = 10**-7;
$PTOL = 10**-4;
$REF_FILE = "output.ref";
$TEST_FILE = "output.dat";
$NDOF = 3;
$RESULT = "cc19.test";

system ("input");
system ("psi3");

open (RE, ">$RESULT") || die "cannot open $RESULT: $!";
select (RE);
printf "\nCC19:\n";

if (abs(seek_nuc($REF_FILE) - seek_nuc($TEST_FILE)) > $TOL) {
  fail_test("Nuclear Repulsion Energy");
}
else {
  pass_test("Nuclear Repulsion Energy");
}

if (abs(seek_scf($REF_FILE) - seek_scf($TEST_FILE)) > $TOL) {
  fail_test("RHF Energy");
}
else {
  pass_test("RHF Energy");
}

if (abs(seek_ccsd($REF_FILE) - seek_ccsd($TEST_FILE)) > $TOL) {
  fail_test("CCSD Energy");
}
else {
  pass_test("CCSD Energy");
}

if (abs(seek_lambda($REF_FILE) - seek_lambda($TEST_FILE)) > $TOL) {
  fail_test("CCSD Lambda Overlap");
}
else {
  pass_test("CCSD Lambda Overlap");
}

@polar_ref = seek_ccsd_polar($REF_FILE);
@polar_test = seek_ccsd_polar($TEST_FILE);

if(!compare_arrays(\@polar_ref,\@polar_test,3,3,$PTOL)) {
  fail_test("CCSD Polarizability");
}
else {
  pass_test("CCSD Polarizability");
}
close(RE);
