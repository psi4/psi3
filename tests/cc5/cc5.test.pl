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
  $REF_FILE = "$SRC_PATH/output.ref";
}
else {
  $REF_FILE = "output.ref";
}
$TEST_FILE = "output.dat";
$RESULT = "cc5.test";

system ("$PSICMD");

$FAIL = 0;

open(RE, ">$RESULT") || die "cannot open $RESULT: $!";
select (RE);
printf "CC5:\n";

if (abs(seek_nuc($REF_FILE) - seek_nuc($TEST_FILE)) > $TOL) {
  fail_test("Nuclear Repulsion Energy"); $FAIL = 1;
}
else {
  pass_test("Nuclear Repulsion Energy");
}

if (abs(seek_scf($REF_FILE) - seek_scf($TEST_FILE)) > $TOL) {
  fail_test("SCF Energy"); $FAIL = 1;
}
else { 
  pass_test("SCF Energy");
}

if (abs(seek_ccsd($REF_FILE) - seek_ccsd($TEST_FILE)) > $TOL) {
  fail_test("CCSD Energy"); $FAIL = 1;
}
else { 
  pass_test("CCSD Energy");
}

if (abs(seek_ccsd_t($REF_FILE) - seek_ccsd_t($TEST_FILE)) > $TOL) {
  fail_test("CCSD(T) Energy"); $FAIL = 1;
}
else {
  pass_test("CCSD(T) Energy");
}

close (RE);

system("cat $RESULT");

exit($FAIL);

