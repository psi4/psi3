#!/usr/bin/perl  

use Cwd;

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

# Test for psi3 -f input -o output -p prefix  (test1)
# and      psi3 input output prefix           (test2)

$OUT_FILE1= "test1.out";
$OUT_FILE2= "test2.out";

$CHK_FILE1= "test1.32";
$CHK_FILE2= "test2.32";

if($SRC_PATH ne "") {
  $IN_FILE = "$SRC_PATH/input.dat";
}

$RESULT = "psi_start.test";

$PSICMD1 = "";
$PSICMD2 = "";

if($EXEC_PATH ne "") {
  $PSICMD1 = "PATH=$EXEC_PATH:\$PATH;export PATH;psi3 -f $IN_FILE -o $OUT_FILE1 -p test1";
  $PSICMD2 = "PATH=$EXEC_PATH:\$PATH;export PATH;psi3 $IN_FILE $OUT_FILE2 test2";
}
else {
  $PSICMD1 = "psi3 -f $IN_FILE -o $OUT_FILE1 -p test1";
  $PSICMD2 = "psi3 $IN_FILE $OUT_FILE2 test2";
}

if($QUIET == 1) {
  $PSICMD1 = $PSICMD1 . " 1>/dev/null 2>/dev/null";
  $PSICMD2 = $PSICMD2 . " 1>/dev/null 2>/dev/null";
}

system ("$PSICMD1");
system ("$PSICMD2");

$testout1 = 0;
$testout2 = 0;
$testchk1 = 0;
$testchk2 = 0;
$FAIL = 0;

$dir = cwd();
opendir(DIR, $dir) || die "cannot open $dir: $!";
@files = readdir(DIR);
for $files (@files) {
  if($files eq "$OUT_FILE1") { $testout1 = 1; }
  if($files eq "$OUT_FILE2") { $testout2 = 1; }
  if($files eq "$CHK_FILE1") { $testchk1 = 1; }
  if($files eq "$CHK_FILE2") { $testchk2 = 1; }
} 
close (DIR);

open(RE, ">$RESULT") || die "cannot open $RESULT: $!";
select (RE);
printf "PSI-START:\n";

if($testout1 == 1 && $testchk1 == 1) { 
  printf "CHECKING psi3 -f input -o output -p prefix     ...PASSED\n";
}
else { 
  printf "CHECKING psi3 -f input -o output -p prefix     ...FAILED\n";
  $FAIL = 1;
}

if($testout2 == 1 && $testchk2 == 1) { 
  printf "CHECKING psi3 input output prefix              ...PASSED\n";
}
else { 
  printf "CHECKING psi3 input output prefix              ...FAILED\n";
  $FAIL = 1;
}

close (RE);

system("cat $RESULT");

exit($FAIL);

