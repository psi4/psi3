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

$TOL = 10**-7;
$PTOL = 10**-4;
if($SRC_PATH ne "") {
  $REF_FILE = "$SRC_PATH/output.ref";
}
else {
  $REF_FILE = "output.ref";
}
$TEST_FILE = "output.dat";
$NDOF = 3;
$RESULT = "cc18.test";

system ("$PSICMD");

$FAIL = 0;

open (RE, ">$RESULT") || die "cannot open $RESULT: $!";
select (RE);
printf "CC18:\n";

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

if (abs(seek_ccsd($REF_FILE) - seek_ccsd($TEST_FILE)) > $TOL) {
  fail_test("CCSD Energy"); $FAIL = 1;
}
else {
  pass_test("CCSD Energy");
}

if (abs(seek_lambda($REF_FILE) - seek_lambda($TEST_FILE)) > $TOL) {
  fail_test("CCSD Lambda Overlap"); $FAIL = 1;
}
else {
  pass_test("CCSD Lambda Overlap");
}

@polar_ref = seek_ccsd_polar($REF_FILE);
@polar_test = seek_ccsd_polar($TEST_FILE);

if(!compare_arrays(\@polar_ref,\@polar_test,3,3,$PTOL)) {
  fail_test("CCSD Polarizability"); $FAIL = 1;
}
else {
  pass_test("CCSD Polarizability");
}
close(RE);

system("cat $RESULT");

exit($FAIL);

