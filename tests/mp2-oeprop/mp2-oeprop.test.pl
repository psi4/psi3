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
$MTOL = 10**-5;
if($SRC_PATH ne "") {
  $REF_FILE = "$SRC_PATH/output.ref";
}
else {
  $REF_FILE = "output.ref";
}
$TEST_FILE = "output.dat";
$RESULT = "mp2-oeprop.test";

system ("$PSICMD");

$FAIL = 0;

open(RE, ">$RESULT") || die "cannot open $RESULT: $!";
select (RE);
printf "MP2-OEPROP:\n";

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

if (abs(seek_mp2($REF_FILE) - seek_mp2($TEST_FILE)) > $TOL) {
  fail_test("MP2 Energy"); $FAIL = 1;
}
else { 
  pass_test("MP2 Energy");
}

@gop_ref = seek_mulliken_gop($REF_FILE);
@gop_test = seek_mulliken_gop($TEST_FILE);

if (!compare_arrays(\@gop_ref,\@gop_test,20,20,$MTOL)) {
  fail_test("Gross orbital populations"); $FAIL = 1;
}
else { 
  pass_test("Gross orbital populations");
}

@abp_ref = seek_mulliken_abp($REF_FILE);
@abp_test = seek_mulliken_abp($TEST_FILE);

if (!compare_arrays(\@abp_ref,\@abp_test,2,2,$MTOL)) {
  fail_test("Atomic bond populations"); $FAIL = 1;
}
else { 
  pass_test("Atomic bond populations");
}

@apnc_ref = seek_mulliken_apnc($REF_FILE);
@apnc_test = seek_mulliken_apnc($TEST_FILE);

if (!compare_arrays(\@apnc_ref,\@apnc_test,2,2,$MTOL)) {
  fail_test("Gross atomic populations and net charges"); $FAIL = 1;
}
else { 
  pass_test("Gross atomic poplulations and net charges");
}

@dipole_ref = seek_dipole($REF_FILE);
@dipole_test = seek_dipole($TEST_FILE);

if (!compare_arrays(\@dipole_ref,\@dipole_test,4,4,$MTOL)) {
  fail_test("Electric dipole moment"); $FAIL = 1;
}
else { 
  pass_test("Electric dipole moment");
}

@angmom_ref = seek_angmom($REF_FILE);
@angmom_test = seek_angmom($TEST_FILE);

if (!compare_arrays(\@angmom_ref,\@angmom_test,3,3,$MTOL)) {
  fail_test("Electronic angular momentum"); $FAIL = 1;
}
else { 
  pass_test("Electronic angular momentum");
}

@epef_ref = seek_epef($REF_FILE);
@epef_test = seek_epef($TEST_FILE);

if (!compare_arrays(\@epef_ref,\@epef_test,8,8,$MTOL)) {
  fail_test("Electrostatic potential and electric field"); $FAIL = 1;
}
else { 
  pass_test("Electrostatic potential and electric field");
}

@edensity_ref = seek_edensity($REF_FILE);
@edensity_test = seek_edensity($TEST_FILE);

if (!compare_arrays(\@edensity_ref,\@edensity_test,2,2,$MTOL)) {
  fail_test("Electron density"); $FAIL = 1;
}
else { 
  pass_test("Electron density");
}

if (abs(seek_mvd($REF_FILE) - seek_mvd($TEST_FILE)) > $TOL) {
  fail_test("Total one-electron MVD terms"); $FAIL = 1;
}
else { 
  pass_test("Total one-electron MVD terms");
}

close (RE);

system("cat $RESULT");

exit($FAIL);

