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
$RESULT = "dboc-rocisd1.test";

system ($PSICMD);

$FAIL = 0;

open(RE, ">$RESULT") || die "cannot open $RESULT $!";
select (RE);
printf "DBOC-RHF1:\n";

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
    
if (abs(seek_ci($REF_FILE) - seek_ci($TEST_FILE)) > $TOL) {
  fail_test("CISD Energy"); $FAIL = 1;
}
else {
  pass_test("CISD Energy");
}

if (abs(seek_dboc($REF_FILE) - seek_dboc($TEST_FILE)) > $TOL) {
  fail_test("DBOC"); $FAIL = 1;
}
else {
  pass_test("DBOC");
}
    
close(RE);

system("cat $RESULT");

exit($FAIL);