#!/usr/bin/perl  

require("../psitest.pl");

$TOL = 10**-8;
$GTOL = 10**-8;
$REF_FILE = "file11.ref";
$TEST_FILE = "psi.file11.dat";
$RESULT = "cc1.test";

system ("input");
system ("psi3");

$natom = seek_natom_file11($REF_FILE,"CCSD");

open(RE, ">$RESULT") || die "cannot open $RESULT $!"; 
select (RE);
printf "\nCC1:\n";

if(abs(seek_energy_file11($REF_FILE,"CCSD") - seek_energy_file11($TEST_FILE,"CCSD")) > $TOL) {
  fail_test("CCSD energy");
}
else {
  pass_test("CCSD energy");
}

@geom_ref = seek_geom_file11($REF_FILE, "CCSD");
@geom_test = seek_geom_file11($TEST_FILE, "CCSD");
if(!compare_arrays(\@geom_ref, \@geom_test, $natom, 3, $GTOL)) {
  fail_test("CCSD Geometry");
}
else {
  pass_test("CCSD Geometry");
}

@grad_ref = seek_grad_file11($REF_FILE, "CCSD");
@grad_test = seek_grad_file11($TEST_FILE, "CCSD");

if(!compare_arrays(\@grad_ref, \@grad_test, $natom, 3, $GTOL)) {
  fail_test("CCSD Gradient");
}
else {
  pass_test("CCSD Gradient");
}

close (RE);
