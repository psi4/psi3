#!/usr/bin/perl

# Automated Test Case Driver
# Added to CVS 1.10.03
# Added different categories for test cases
# --ci --cc --scf --small --medium --large
# --all --sp --geom --freq --excite
# The default is no argument and all the
# test cases will run except the "--large" ones.
# The arguments --clean and --cleanall
# were also added in case the test cases
# need to be run more than once
# 1.14.03

use Cwd;

if ($ARGV[0] eq "--ci") {
  @DIR_NAMES = ( "cisd-sp",
                 "cisd-opt-numer",
                 "casscf-sp",
                 "cis-sp"
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
                 "scf-opt2",
                 "scf-opt3",
                 "scf-opt4",
                 "scf-opt5",
                 "scf-opt6",
                 "scf-opt7",
                 "scf-opt8",
                 "scf-opt9",
                 "scf-opt10",
                 "scf-opt11",
                 "scf-opt-numer",
                 "scf-freq",
                 "scf-freq-numer",
                 "scf-freq-symm-numer",
                 "scf-polar",
                 "rhf-stab",
                 "uhf-stab",
                 "rohf-stab",
                 "extrema-zmat",
                 "extrema-deloc"
               );
}
elsif ($ARGV[0] eq "--dboc") {
  @DIR_NAMES = ( "dboc-rhf1",
                 "dboc-rohf1",
                 "dboc-uhf1",
                 "dboc-rcisd1",
                 "dboc-rocisd1"
               );
}
elsif ($ARGV[0] eq "--small") {
  @DIR_NAMES = ( "scf-opt",
                 "scf-opt2",
                 "scf-opt3",
                 "scf-opt4",
                 "scf-opt5",
                 "scf-opt6",
                 "scf-opt7",
                 "scf-opt8",
                 "scf-opt9",
                 "scf-opt10",
                 "scf-opt11",
                 "scf-opt-numer",
                 "scf-freq",
                 "scf-freq-numer",
                 "scf-freq-symm-numer",
                 "scf-polar",
                 "rhf-stab",
                 "uhf-stab",
                 "rohf-stab",
		 "mp2-sp",
		 "mp2-direct-sp",
                 "extrema-zmat",
                 "extrema-deloc",
                 "cis-sp",
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
                 "cc19",
                 "dboc-rhf1",
                 "dboc-rohf1",
                 "dboc-uhf1",
                 "mp2r12-sp1"
               );
}
elsif ($ARGV[0] eq "--medium") {
  @DIR_NAMES = ( "cc4",
                 "dboc-rcisd1",
                 "dboc-rocisd1" 
               );
}
elsif ($ARGV[0] eq "--large") {
  @DIR_NAMES = ( "cc5",
                 "cc6",
                 "cc7"
               );
}
elsif ($ARGV[0] eq "--all") {
  @DIR_NAMES = ( "scf-opt",
                 "scf-opt2",
                 "scf-opt3",
                 "scf-opt4",
                 "scf-opt5",
                 "scf-opt6",
                 "scf-opt7",
                 "scf-opt8",
                 "scf-opt9",
                 "scf-opt10",
                 "scf-opt11",
                 "scf-opt-numer",
                 "scf-freq",
                 "scf-freq-numer",
                 "scf-freq-symm-numer",
                 "scf-polar",
                 "rhf-stab",
                 "uhf-stab",
                 "rohf-stab",
		 "mp2-sp",
		 "mp2-direct-sp",
                 "extrema-zmat",
                 "extrema-deloc",
                 "cis-sp",
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
                 "cc19",
                 "dboc-rhf1",
                 "dboc-rohf1",
                 "dboc-uhf1",
                 "dboc-rcisd1",
                 "dboc-rocisd1",
                 "mp2r12-sp1"
               );
}
elsif ($ARGV[0] eq "--sp") {
  @DIR_NAMES = ( "mp2-sp",
                 "mp2-direct-sp",
                 "cisd-sp",
                 "casscf-sp",
                 "cis-sp",
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
                 "cc19",
                 "dboc-rhf1",
                 "dboc-rohf1",
                 "dboc-uhf1",
                 "dboc-rcisd1",
                 "dboc-rocisd1",
                 "mp2r12-sp1"
               );
}
elsif ($ARGV[0] eq "--geom") {
  @DIR_NAMES = ( "scf-opt",
                 "scf-opt2",
                 "scf-opt3",
                 "scf-opt4",
                 "scf-opt5",
                 "scf-opt6",
                 "scf-opt7",
                 "scf-opt8",
                 "scf-opt9",
                 "scf-opt10",
                 "scf-opt11",
                 "scf-opt-numer",
                 "cisd-opt-numer",
                 "cc1",
                 "cc2",
                 "cc13",
                 "cc14",
                 "extrema-zmat",
                 "extrema-deloc"
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
                 "scf-opt2",
                 "scf-opt3",
                 "scf-opt4",
                 "scf-opt5",
                 "scf-opt6",
                 "scf-opt7",
                 "scf-opt8",
                 "scf-opt9",
                 "scf-opt10",
                 "scf-opt11",
                 "scf-opt-numer",
                 "scf-freq",
                 "scf-freq-numer",
                 "scf-freq-symm-numer",
                 "scf-polar",
                 "rhf-stab",
                 "uhf-stab",
                 "rohf-stab",
		 "mp2-sp",
		 "mp2-oeprop",
		 "mp2-direct-sp",
                 "extrema-zmat",
                 "extrema-deloc",
                 "cis-sp",
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
                 "cc19",
                 "dboc-rhf1",
                 "dboc-rohf1",
                 "dboc-uhf1",
                 "dboc-rcisd1",
                 "dboc-rocisd1",
                 "mp2r12-sp1"
#                 "Basis_Tests"
                );
}

$pwd = cwd();

if ($ARGV[0] != "--clean" || $ARGV[1] != "--clean") {
  unlink <test-case-results>;
}

foreach $name (@DIR_NAMES) {
  print "$name\n";
  chdir ("$name");

  if ($ARGV[0] eq "--clean" || $ARGV[1] eq "--clean") {
    unlink <core>;
    unlink <output.dat>;
    unlink <file*.dat>;
    unlink <psi.*>;
    unlink <*.test>;
    unlink <timer.dat>;
    unlink <dboc.findif.out>;
    chdir ("$pwd");
  }
  else {
    system ("./$name.test.pl");
    system ("psiclean");
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
}

