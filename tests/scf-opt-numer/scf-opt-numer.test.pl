#!/usr/bin/perl  

require("../psitest.pl");

$TOL = 10**-8;
$GTOL = 10**-8;
$REF_FILE = "file11.ref";
$TEST_FILE = "psi.file11.dat";
$RESULT = "scf-opt-numer.test";

system ("input");
system ("psi3");

$natom = seek_natom_file11($REF_FILE,"iteration");

open(RE, ">$RESULT") || die "cannot open $RESULT $!"; 
select (RE);
printf "\nSCF-OPT-NUMER:\n";

if(abs(seek_energy_file11($REF_FILE,"iteration") - seek_energy_file11($TEST_FILE,"iteration")) > $TOL) {
  fail_test("SCF energy");
}
else {
  pass_test("SCF energy");
}

@geom_ref = seek_geom_file11($REF_FILE, "iteration");
@geom_test = seek_geom_file11($TEST_FILE, "iteration");
if(!compare_arrays(\@geom_ref, \@geom_test, $natom, 3, $GTOL)) {
  fail_test("SCF Geometry");
}
else {
  pass_test("SCF Geometry");
}

@grad_ref = seek_grad_file11($REF_FILE, "iteration");
@grad_test = seek_grad_file11($TEST_FILE, "iteration");

if(!compare_arrays(\@grad_ref, \@grad_test, $natom, 3, $GTOL)) {
  fail_test("SCF Gradient");
}
else {
  pass_test("SCF Gradient");
}

close (RE);
