#!/usr/bin/perl  

require("../psitest.pl");

$TOL = 10**-10;
$PTOL = 10**-4;
$REF_FILE = "output.ref";
$TEST_FILE = "output.dat";
$NDOF = 3;
$RESULT = "scf-polar.test";

system ("psi3");

open (RE, ">$RESULT") || die "cannot open $RESULT: $!";
select (RE);
printf "\nSCF-POLAR:\n";

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

@polar_ref = seek_scf_polar($REF_FILE);
@polar_test = seek_scf_polar($TEST_FILE);

if(!compare_arrays(\@polar_ref,\@polar_test,3,3,$PTOL)) {
  fail_test("SCF Polarizability");
}
else {
  pass_test("SCF Polarizability");
}
close(RE);
