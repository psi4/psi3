#!/usr/bin/perl  

require ("../psitest.pl");

$TOL = 10**-8;
$REF_FILE = "output.ref";
$TEST_FILE = "output.dat";
$RESULT = "cisd-sp.test";

system ("input");
system ("psi3");

open(RE, ">$RESULT") || die "cannot open $RESULT $!";
select (RE);
printf "\nCISD-SP:\n";

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
    
if (abs(seek_ci($REF_FILE) - seek_ci($TEST_FILE)) > $TOL) {
  fail_test("CISD Energy");
}
else {
  pass_test("CISD Energy");
}
    
close(RE);

