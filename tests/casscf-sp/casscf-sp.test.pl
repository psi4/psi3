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
if($SRC_PATH ne "") {
  $REF_FILE = "$SRC_PATH/file14.ref";
}
else {
  $REF_FILE = "file14.ref";
}
$TEST_FILE = "file14.dat";
if($SRC_PATH ne "") {
  $REF_OUT = "$SRC_PATH/output.ref";
}
else {
  $REF_OUT = "output.ref";
}
$TEST_OUT = "output.dat";
$RESULT = "casscf-sp.test";

system ($PSICMD);

$FAIL = 0;

open(RE, ">$RESULT") || die "cannot open $RESULT $!";
select (RE);
printf "CASSCF-SP:\n";

if (abs(seek_nuc($REF_OUT) - seek_nuc($TEST_OUT)) > $TOL) {
  fail_test("Nuclear Repulsion Energy"); $FAIL = 1;
}
else {
  pass_test("Nuclear Repulsion Energy");
}

if (abs(seek_scf($REF_OUT) - seek_scf($TEST_OUT)) > $TOL) {
  fail_test("RHF Energy"); $FAIL = 1;
}
else {
  pass_test("RHF Energy");
}
    
if (abs(seek_casscf($REF_FILE) - seek_casscf($TEST_FILE)) > $TOL) {
  fail_test("CASSCF Energy"); $FAIL = 1;
}
else {
  pass_test("CASSCF Energy");
}
    
close(RE);

system("cat $RESULT");

exit($FAIL);
