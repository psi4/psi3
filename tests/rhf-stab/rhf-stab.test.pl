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
$STAB_TOL = 10**-4;
if($SRC_PATH ne "") {
  $REF_FILE = "$SRC_PATH/output.ref";
}
else {
  $REF_FILE = "output.ref";
}
$TEST_FILE = "output.dat";
$RESULT = "rhf-stab.test";
$NUM_EVALS = 5;
$NUM_SYMMS = 4;

system ($PSICMD);

$FAIL = 0;

open(RE, ">$RESULT") || die "cannot open $RESULT: $!";
select (RE);
printf "RHF-STAB:\n";

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

@stab_ref = seek_stab($REF_FILE, "RHF->RHF", $NUM_EVALS, $NUM_SYMMS);
@stab_test = seek_stab($TEST_FILE, "RHF->RHF", $NUM_EVALS, $NUM_SYMMS);

if(!compare_arrays(\@stab_ref,\@stab_test,$NUM_EVALS,$NUM_SYMMS,$STAB_TOL))
{
  fail_test("RHF->RHF stability"); $FAIL = 1;
}
else {
  pass_test("RHF->RHF stability");
}

@stab_ref = seek_stab($REF_FILE, "RHF->UHF", $NUM_EVALS, $NUM_SYMMS);
@stab_test = seek_stab($TEST_FILE, "RHF->UHF", $NUM_EVALS, $NUM_SYMMS);

if(!compare_arrays(\@stab_ref,\@stab_test,$NUM_EVALS,$NUM_SYMMS,$STAB_TOL))
{
  fail_test("RHF->UHF stability"); $FAIL = 1;
}
else {
  pass_test("RHF->UHF stability");
}

close (RE);

system("cat $RESULT");

exit($FAIL);
