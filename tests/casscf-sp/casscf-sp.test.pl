#!/usr/bin/perl  

require ("../psitest.pl");

$TOL = 10**-8;
$REF_FILE = "file14.ref";
$TEST_FILE = "file14.dat";
$REF_OUT = "output.ref";
$TEST_OUT = "output.dat";
$RESULT = "casscf-sp.test";

system ("psi3");

open(RE, ">$RESULT") || die "cannot open $RESULT $!";
select (RE);
printf "\nCASSCF-SP:\n";

if (abs(seek_nuc($REF_OUT) - seek_nuc($TEST_OUT)) > $TOL) {
  fail_test("Nuclear Repulsion Energy");
}
else {
  pass_test("Nuclear Repulsion Energy");
}

if (abs(seek_scf($REF_OUT) - seek_scf($TEST_OUT)) > $TOL) {
  fail_test("RHF Energy");
}
else {
  pass_test("RHF Energy");
}
    
if (abs(seek_casscf($REF_FILE) - seek_casscf($TEST_FILE)) > $TOL) {
  fail_test("CASSCF Energy");
}
else {
  pass_test("CASSCF Energy");
}
    
close(RE);

