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
$RESULT = "scf-freq-first-5.test";
$NDOF = 3;

system ($PSICMD);

$FAIL = 0;

$natom = seek_natom_file11($REF_FILE,"SCF");

open(RE, ">$RESULT") || die "cannot open $RESULT $!"; 
select (RE);
printf "SCF-FREQ-NUMER:\n";

if(abs(seek_energy_file11($REF_FILE,"SCF") - seek_energy_file11($TEST_FILE,"SCF")) > $TOL) {
  fail_test("SCF energy"); $FAIL = 1;
}
else {
  pass_test("SCF energy");
}

@geom_ref = seek_geom_file11($REF_FILE, "SCF");
@geom_test = seek_geom_file11($TEST_FILE, "SCF");
if(!compare_arrays(\@geom_ref, \@geom_test, $natom, 3, $GTOL)) {
  fail_test("SCF Geometry"); $FAIL = 1;
}
else {
  pass_test("SCF Geometry");
}

@grad_ref = seek_grad_file11($REF_FILE, "SCF");
@grad_test = seek_grad_file11($TEST_FILE, "SCF");

if(!compare_arrays(\@grad_ref, \@grad_test, $natom, 3, $GTOL)) {
  fail_test("SCF Gradient"); $FAIL = 1;
}
else {
  pass_test("SCF Gradient");
}

@freq_ref = seek_findif_freq($REF_OUT,"Harmonic Vibrational Frequencies",$NDOF);
@freq_test = seek_findif_freq($TEST_OUT,"Harmonic Vibrational Frequencies",$NDOF);

if(!compare_arrays(\@freq_ref, \@freq_test, $NDOF, 1, $HTOL)) {
  fail_test("SCF Frequencies"); $FAIL = 1;
}
else {
  pass_test("SCF Frequencies");
}    

close (RE);

system("cat $RESULT");

exit($FAIL);
