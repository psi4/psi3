#!/usr/bin/perl  

while ($ARGV = shift) {
  if   ("$ARGV" eq "-q") { $QUIET = 1; }
  elsif("$ARGV" eq "-i") { $SRC_PATH = shift; }
  elsif("$ARGV" eq "-x") { $EXEC_PATH = shift; }
}

if($SRC_PATH ne "") {
  require($SRC_PATH . "/../psitest.pl");
}
else {
  require("../psitest.pl");
}

# build the command for the psi3 driver
$PSICMD = build_psi_cmd($QUIET, $SRC_PATH, $EXEC_PATH);

$TOL = 10**-8;
$GTOL = 10**-8;
$HTOL = 10**-2;
if($SRC_PATH ne "") {
  $REF_FILE = "$SRC_PATH/file11.ref";
}
else {
  $REF_FILE = "file11.ref";
}
$TEST_FILE = "psi.file11.dat";
if($SRC_PATH ne "") {
  $REF_OUT = "$SRC_PATH/output.ref";
}
else {
  $REF_OUT = "output.ref";
}
$TEST_OUT = "output.dat";
$RESULT = "cc3.test";
$NDOF = 3;

system ("$PSICMD");

$FAIL = 0;

$natom = seek_natom_file11($REF_FILE,"CCSD");

open(RE, ">$RESULT") || die "cannot open $RESULT $!";
select (RE);
printf "CC3:\n";

if(abs(seek_energy_file11($REF_FILE,"CCSD") - seek_energy_file11($TEST_FILE,"CCSD")) > $TOL) {
  fail_test("CCSD energy"); $FAIL = 1;
}
else {
  pass_test("CCSD energy");
}

@geom_ref = seek_geom_file11($REF_FILE, "CCSD");
@geom_test = seek_geom_file11($TEST_FILE, "CCSD");
if(!compare_arrays(\@geom_ref, \@geom_test, $natom, 3, $GTOL)) {
  fail_test("CCSD Geometry"); $FAIL = 1;
}
else {
  pass_test("CCSD Geometry");
}

@grad_ref = seek_grad_file11($REF_FILE, "CCSD");
@grad_test = seek_grad_file11($TEST_FILE, "CCSD");

if(!compare_arrays(\@grad_ref, \@grad_test, $natom, 3, $GTOL)) {
  fail_test("CCSD Gradient"); $FAIL = 1;
}
else {
  pass_test("CCSD Gradient");
}

@freq_ref = seek_findif_freq($REF_OUT,"Harmonic Vibrational Frequencies",$NDOF);
@freq_test = seek_findif_freq($TEST_OUT,"Harmonic Vibrational Frequencies",$NDOF);

if(!compare_arrays(\@freq_ref, \@freq_test, $NDOF, 1, $HTOL)) {
  fail_test("CCSD Frequencies"); $FAIL = 1;
}
else {
  pass_test("CCSD Frequencies");
}

close (RE);

system("cat $RESULT");

exit($FAIL);

