#!/usr/bin/perl

use Cwd;

@DIR_NAMES = ( "scf-opt",
               "scf-opt-numer",
#               "scf-freq",
#               "scf-freq-numer",
#               "scf-freq-symm-numer",
#               "cisd-sp",
#               "cisd-opt-numer",
#               "casscf-sp",
#               "cc1", 
#               "cc2",
#               "cc3",
#               "cc4",
#               "cc5",
#               "cc6",
#               "cc7",
#               "cc8",
#               "cc9",
#               "cc10",
#               "cc11",
#               "cc12",
#               "Basis_Tests"
                );

$pwd = cwd();

foreach $name (@DIR_NAMES) {
  chdir ("$name");
  system ("./$name.test.pl");
  open (T,"$name.test");
  @data = <T>;
  chdir ("$pwd");
  open (F,">>test-case-results");
  select (F);
  seek(Q,0,2);
  print @data;
  close (T);
  close (F);
}

