#!/usr/bin/perl  

require("../psitest.pl");

$TOL = 10**-8;
$REF_FILE = "output.ref";
$TEST_FILE = "output.dat";
$RESULT = "cc8.test";

system ("input");
system ("psi3");

open(RE, ">$RESULT") || die "cannot open $RESULT: $!";
select (RE);
printf "\nCC8:\n";

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

if (abs(seek_ccsd_t($REF_FILE) - seek_ccsd_t($TEST_FILE)) > $TOL) {
  fail_test("CCSD(T) Energy");
}
else {
  pass_test("CCSD(T) Energy");
}

close (RE);
