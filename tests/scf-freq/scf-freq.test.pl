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
$HTOL = 10**-2;
$ITOL = 10**-3;
$NDOF = 9;
if($SRC_PATH ne "") {
  $REF_FILE = "$SRC_PATH/output.ref";
}
else {
  $REF_FILE = "output.ref";
}
$TEST_FILE = "output.dat";
$RESULT = "scf-freq.test";

system ($PSICMD);

$FAIL = 0;

open(RE, ">$RESULT") || die "cannot open $RESULT: $!";
select (RE);
printf "SCF-FREQ:\n";

if (abs(seek_nuc($REF_FILE) - seek_nuc($TEST_FILE)) > $TOL) {
  fail_test("Nuclear Repulsion Energy"); $FAIL = 1;
}
else {
  pass_test("Nuclear Repulsion Energy");
}

if (abs(seek_scf($REF_FILE) - seek_scf($TEST_FILE)) > $TOL) {
  fail_test("RHF Energy"); $FAIL = 1;
}
else {
  pass_test("RHF Energy");
}

@freq_ref = seek_anal_freq($REF_FILE,"Harmonic Frequency",$NDOF);
@freq_test = seek_anal_freq($TEST_FILE,"Harmonic Frequency",$NDOF);

if(!compare_arrays(\@freq_ref, \@freq_test, $NDOF, 1, $HTOL)) {
  fail_test("SCF Frequencies"); $FAIL = 1;
}
else {
  pass_test("SCF Frequencies");
}

@int_ref = seek_int($REF_FILE,"Harmonic Frequency",$NDOF);
@int_test = seek_int($TEST_FILE,"Harmonic Frequency",$NDOF);

if(!compare_arrays(\@int_ref, \@int_test, $NDOF, 1, $ITOL)) {
  fail_test("SCF Intensities"); $FAIL = 1;
}
else {
  pass_test("SCF Intensities");
}
    
close(RE);

system("cat $RESULT");

exit($FAIL);
