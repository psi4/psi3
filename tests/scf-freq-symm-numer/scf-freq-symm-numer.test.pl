#!/usr/bin/perl  

require("../psitest.pl");

$TOL = 10**-8;
$GTOL = 10**-8;
$HTOL = 10**-2;
$REF_FILE = "file11.ref";
$TEST_FILE = "psi.file11.dat";
$REF_OUT = "output.ref";
$TEST_OUT = "output.dat";
$RESULT = "scf-freq-symm-numer.test";
$NDOF=2;

system ("psi3");

$natom = seek_natom_file11($REF_FILE,"SCF");

open(RE, ">$RESULT") || die "cannot open $RESULT $!"; 
select (RE);
printf "\nSCF-FREQ-SYMM-NUMER:\n";

if(abs(seek_energy_file11($REF_FILE,"SCF") - seek_energy_file11($TEST_FILE,"SCF")) > $TOL) {
  fail_test("SCF energy");
}
else {
  pass_test("SCF energy");
}

@geom_ref = seek_geom_file11($REF_FILE, "SCF");
@geom_test = seek_geom_file11($TEST_FILE, "SCF");
if(!compare_arrays(\@geom_ref, \@geom_test, $natom, 3, $GTOL)) {
  fail_test("SCF Geometry");
}
else {
  pass_test("SCF Geometry");
}

@grad_ref = seek_grad_file11($REF_FILE, "SCF");
@grad_test = seek_grad_file11($TEST_FILE, "SCF");

if(!compare_arrays(\@grad_ref, \@grad_test, $natom, 3, $GTOL)) {
  fail_test("SCF Gradient");
}
else {
  pass_test("SCF Gradient");
}

@freq_ref = seek_findif_freq($REF_OUT,"Harmonic Vibrational Frequencies",$NDOF);
@freq_test = seek_findif_freq($TEST_OUT,"Harmonic Vibrational Frequencies",$NDOF);

if(!compare_arrays(\@freq_ref, \@freq_test, $NDOF, 1, $HTOL)) {
  fail_test("SCF Frequencies");
}
else {
  pass_test("SCF Frequencies");
}    

close (RE);
