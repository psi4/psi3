#!/usr/bin/perl  

require("../psitest.pl");

$TOL = 10**-8;
$REF_FILE = "output.ref";
$TEST_FILE = "output.dat";
$RESULT = "cc16.test";

system ("input");
system ("psi3");

open(RE, ">$RESULT") || die "cannot open $RESULT $!"; 
select (RE);
printf "\nCC16:\n";

if (abs (seek_scf($REF_FILE) - seek_scf($TEST_FILE)) > $TOL) {
  fail_test("SCF Energy");
}
else {
  pass_test("SCF Energy");
}

if (abs (seek_bccd($REF_FILE) - seek_bccd($TEST_FILE)) > $TOL) {
  fail_test("B-CCD Energy");
}
else {
  pass_test("B-CCD Energy");
}

if (abs (seek_ccsd_t($REF_FILE) - seek_ccsd_t($TEST_FILE)) > $TOL) {
  fail_test("B-CCD(T) Energy");
}
else {
  pass_test("B-CCD(T) Energy");
}
close (RE);
