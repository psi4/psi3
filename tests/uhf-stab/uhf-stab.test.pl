#!/usr/bin/perl  

require("../psitest.pl");

$TOL = 10**-8;
$STAB_TOL = 10**-4;
$REF_FILE = "output.ref";
$TEST_FILE = "output.dat";
$RESULT = "uhf-stab.test";
$NUM_EVALS = 5;
$NUM_SYMMS = 4;

system ("input");
system ("psi3");

open(RE, ">$RESULT") || die "cannot open $RESULT: $!";
select (RE);
printf "\nUHF-STAB:\n";

if (abs(seek_nuc($REF_FILE) - seek_nuc($TEST_FILE)) > $TOL) {
  fail_test("Nuclear Repulsion Energy");
}
else {
  pass_test("Nuclear Repulsion Energy");
}

if (abs(seek_scf($REF_FILE) - seek_scf($TEST_FILE)) > $TOL) {
  fail_test("SCF Energy");
}
else {
  pass_test("SCF Energy");
}

@stab_ref = seek_stab($REF_FILE, "UHF->UHF", $NUM_EVALS, $NUM_SYMMS);
@stab_test = seek_stab($TEST_FILE, "UHF->UHF", $NUM_EVALS, $NUM_SYMMS);

if(!compare_arrays(\@stab_ref,\@stab_test,$NUM_EVALS,$NUM_SYMMS,$STAB_TOL))
{
  fail_test("UHF->UHF stability");
}
else {
  pass_test("UHF->UHF stability");
}

close (RE);
