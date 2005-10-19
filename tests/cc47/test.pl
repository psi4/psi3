#!/usr/bin/perl  

require("../psitest.pl");

$PSICMD = "psi3";

$TOL = 10**-8;
$TTOL = 10**-4;

$REF_FILE = "output.ref";
$TEST_FILE = "output.dat";
$RESULT = "output.test";

system ($PSICMD);

$FAIL = 0;

open(RE, ">$RESULT") || die "cannot open $RESULT: $!";
select (RE);
printf "EOM CCSD OEPROP:\n";

$NSTATES = seek_nstates($REF_FILE);

$LABEL = "Nuclear Repulsion Energy";
if (abs(seek_nuc($REF_FILE) - seek_nuc($TEST_FILE)) > $TOL) {
  printf "%-70s...FAILED\n",$LABEL;
  $FAIL = 1;
}
else {
  printf "%-70s...PASSED\n",$LABEL;
}

$LABEL = "RHF Energy";
if (abs(seek_scf($REF_FILE) - seek_scf($TEST_FILE)) > $TOL) {
  printf "%-70s...FAILED\n",$LABEL;
  $FAIL = 1;
}
else {
  printf "%-70s...PASSED\n",$LABEL;
}

@int_ref = seek_excitation_energy($REF_FILE,"OS       RS",$NSTATES);
@int_test = seek_excitation_energy($TEST_FILE,"OS       RS",$NSTATES);

$LABEL = "Excitation Energy";
if(!compare_arrays(\@int_ref, \@int_test, $NSTATES, 1, $TOL)) {
  printf "%-70s...FAILED\n",$LABEL;
  $FAIL = 1;
}
else {
  printf "%-70s...PASSED\n",$LABEL;
}
    
@int_ref = seek_osc_str($REF_FILE,"OS       RS",$NSTATES);
@int_test = seek_osc_str($TEST_FILE,"OS       RS",$NSTATES);

$LABEL = "Oscillator Strength";
if(!compare_arrays(\@int_ref, \@int_test, $NSTATES, 1, $TTOL)) {
  printf "%-70s...FAILED\n",$LABEL;
  $FAIL = 1;
}
else {
  printf "%-70s...PASSED\n",$LABEL;
}
    
@int_ref = seek_rot_str($REF_FILE,"OS       RS",$NSTATES);
@int_test = seek_rot_str($TEST_FILE,"OS       RS",$NSTATES);

$LABEL = "Rotational Strength";
if(!compare_arrays(\@int_ref, \@int_test, $NSTATES, 1, $TTOL)) {
  printf "%-70s...FAILED\n",$LABEL;
  $FAIL = 1;
}
else {
  printf "%-70s...PASSED\n",$LABEL;
}
    
close(RE);

system("cat $RESULT");

exit($FAIL);

