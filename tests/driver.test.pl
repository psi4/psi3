#!/usr/bin/perl

# Automated Test Case Driver
# Added to CVS 1.10.03
# Added different categories for test cases
# --ci --cc --scf --small --medium --large
# --standard --sp --geom --freq --excite
# the default is no argument and all the
# test cases will run
# The arguments --clean and --cleanall
# were also added in case the test cases
# need to be run more than once
# 1.14.03

use Cwd;

if ($ARGV[0] eq "--ci") {
  @DIR_NAMES = ( "cisd-sp",
                 "cisd-opt-numer",
                 "casscf-sp"
               );
}
elsif ($ARGV[0] eq "--cc") {
  @DIR_NAMES = ( "cc1",
                 "cc2",
                 "cc3",
                 "cc4",
                 "cc5",
                 "cc6",
                 "cc7",
                 "cc8",
                 "cc9",
                 "cc10",
                 "cc11",
                 "cc12",
                 "cc13",
                 "cc14",
                 "cc15",
                 "cc16",
                 "cc17",
                 "cc18",
                 "cc19"
               );
}
elsif ($ARGV[0] eq "--scf") {
  @DIR_NAMES = ( "scf-opt",
                 "scf-opt-numer",
                 "scf-freq",
                 "scf-freq-numer",
                 "scf-freq-symm-numer",
                 "scf-polar",
                 "rhf-stab",
                 "uhf-stab",
                 "rohf-stab"
               );
}
elsif ($ARGV[0] eq "--small") {
  @DIR_NAMES = ( "scf-opt",
                 "scf-opt-numer",
                 "scf-freq",
                 "scf-freq-numer",
                 "scf-freq-symm-numer",
                 "scf-polar",
                 "rhf-stab",
                 "uhf-stab",
                 "rohf-stab",
                 "cisd-sp",
                 "cisd-opt-numer",
                 "casscf-sp",
                 "cc1",
                 "cc2",
                 "cc3",
                 "cc8",
                 "cc9",
                 "cc10",
                 "cc11",
                 "cc12",
                 "cc13",
                 "cc14",
                 "cc15",
                 "cc16",
                 "cc17",
                 "cc18",
                 "cc19"
               );
}
elsif ($ARGV[0] eq "--medium") {
  @DIR_NAMES = ( "cc4" 
               );
}
elsif ($ARGV[0] eq "--large") {
  @DIR_NAMES = ( "cc5",
                 "cc6",
                 "cc7"
               );
}
elsif ($ARGV[0] eq "--standard") {
  @DIR_NAMES = ( "scf-opt",
                 "scf-opt-numer",
                 "scf-freq",
                 "scf-freq-numer",
                 "scf-freq-symm-numer",
                 "scf-polar",
                 "rhf-stab",
                 "uhf-stab",
                 "rohf-stab",
                 "cisd-sp",
                 "cisd-opt-numer",
                 "casscf-sp",
                 "cc1",
                 "cc2",
                 "cc3",
                 "cc4",
                 "cc8",
                 "cc9",
                 "cc10",
                 "cc11",
                 "cc12",
                 "cc13",
                 "cc14",
                 "cc15",
                 "cc16",
                 "cc17",
                 "cc18",
                 "cc19"
               );
}
elsif ($ARGV[0] eq "--sp") {
  @DIR_NAMES = ( "cisd-sp",
                 "casscf-sp",
                 "cc4",
                 "cc5",
                 "cc6",
                 "cc7",
                 "cc8",
                 "cc9",
                 "cc10",
                 "cc11",
                 "cc15",
                 "cc16",
                 "cc17",
                 "cc18",
                 "cc19"
               );
}
elsif ($ARGV[0] eq "--geom") {
  @DIR_NAMES = ( "scf-opt",
                 "scf-opt-numer",
                 "cisd-opt-numer",
                 "cc1",
                 "cc2",
                 "cc13",
                 "cc14"
               );
}
elsif ($ARGV[0] eq "--freq") {
  @DIR_NAMES = ( "scf-freq",
                 "scf-freq-numer",
                 "scf-freq-symm-numer",
                 "cc3"
               );
}
elsif ($ARGV[0] eq "--excite") {
  @DIR_NAMES = ( "cc12", "cc17" );
}
else {
  @DIR_NAMES = ( "scf-opt",
                 "scf-opt-numer",
                 "scf-freq",
                 "scf-freq-numer",
                 "scf-freq-symm-numer",
                 "scf-polar",
                 "rhf-stab",
                 "uhf-stab",
                 "rohf-stab",
                 "cisd-sp",
                 "cisd-opt-numer",
                 "casscf-sp",
                 "cc1", 
                 "cc2",
                 "cc3",
                 "cc4",
                 "cc5",
                 "cc6",
                 "cc7",
                 "cc8",
                 "cc9",
                 "cc10",
                 "cc11",
                 "cc12",
                 "cc13",
                 "cc14",
                 "cc15",
                 "cc16",
                 "cc17",
                 "cc18",
                 "cc19"
#                 "Basis_Tests"
                );
}

$pwd = cwd();

foreach $name (@DIR_NAMES) {
  chdir ("$name");
  system ("./$name.test.pl");
  system ("psiclean");
  open (T,"$name.test");
  @data = <T>;

  if ($ARGV[0] eq "--clean" || $ARGV[1] eq "--clean") {
    unlink <psi.*>;
    unlink <*.test>;
    unlink <timer.dat>;
  }
  elsif ($ARGV[0] eq "--cleanall" || $ARGV[1] eq "--cleanall") {
    unlink <output.dat>;
    unlink <file*.dat>;
    unlink <psi.*>;
    unlink <*.test>;
    unlink <timer.dat>;
  }

  chdir ("$pwd");
  open (F,">>test-case-results");
  select (F);
  seek(Q,0,2);
  print @data;
  close (T);
  close (F);
}

