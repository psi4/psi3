#!/usr/bin/perl  

require("../psitest.pl");

$TOL = 10**-8;
$GTOL = 10**-8;
$REF_FILE = "file11.ref";
$TEST_FILE = "psi.file11.dat";
$RESULT = "cisd-opt-numer.test";

system ("psi3");

$natom = seek_natom_file11($REF_FILE,"iteration");

open(RE, ">$RESULT") || die "cannot open $RESULT $!"; 
select (RE);
printf "\nCISD-OPT-NUMER:\n";

if(abs(seek_energy_file11($REF_FILE,"interation") - seek_energy_file11($TEST_FILE,"iteration")) > $TOL) {
  fail_test("CISD energy");
}
else {
  pass_test("CISD energy");
}

@geom_ref = seek_geom_file11($REF_FILE, "iteration");
@geom_test = seek_geom_file11($TEST_FILE, "iteration");
if(!compare_arrays(\@geom_ref, \@geom_test, $natom, 3, $GTOL)) {
  fail_test("CISD Geometry");
}
else {
  pass_test("CISD Geometry");
}

@grad_ref = seek_grad_file11($REF_FILE, "iteration");
@grad_test = seek_grad_file11($TEST_FILE, "iteration");

if(!compare_arrays(\@grad_ref, \@grad_test, $natom, 3, $GTOL)) {
  fail_test("CISD Gradient");
}
else {
  pass_test("CISD Gradient");
}

close (RE);
