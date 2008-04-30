#!/usr/bin/perl

#
# Prepare PSI code for shipment to non-developer parties
# Edward Valeev (evaleev at vt.edu) 2008
#

my $svnurl = "https://sirius.chem.vt.edu/svn/psi3/trunk";

my $nargs = $#ARGV;
usage() if ($nargs != 0);

my $rev = shift;
# if asked for HEAD, determine the actual revision number from the repository
if ($rev eq "HEAD") {
  $_ = `svn info $svnurl | grep 'Revision:'`;
  chomp; my @tokens = split;
  $rev = $tokens[1];
}

my $briefdirname = "psi3-r$rev";
system("rm -rf $briefdirname");
system("svn export -r $rev $svnurl $briefdirname");

# extract psi version and buildid from checked-out configure.ac
my $version;  my $buildid;
{
  $_ = `grep psi_version $briefdirname/configure.ac | grep define`;
  chomp;  s/\[/ /g;  s/\]/ /g;  my @tokens = split; $version = $tokens[3];
  $_ = `grep psi_buildid $briefdirname/configure.ac | grep define`;
  chomp;  s/\[/ /g;  s/\]/ /g;  @tokens = split; $buildid = $tokens[3];
}

# form the full name and rename the directory
my $dirname = "psi-$version-$buildid-r$rev";
system("mv $briefdirname $dirname");

# sanitize
sanitize();

# prepare for build
system ("(cd $dirname; autoconf; rm -rf autom4te.cache)");

# make .tgz
system("rm -rf $dirname.tar.gz");
system("tar cvzf $dirname.tar.gz $dirname");
system("rm -rf $dirname");

print "Prepared $dirname.tar.gz: psi3 version $version buildid $buildid revision $rev\n";

exit;

0;

sub usage() {

  printf STDOUT "Usage: export_trunk.pl REV\n";
  printf STDOUT "where: REV -- the revision to export (numeric revision specifier or HEAD)\n";

}

#
# Rip out pieces not ready for distribution
#
sub sanitize() {

  system("rm -rf $dirname/misc");
  
  my @binaries = ("psimrcc", "mcscf", "nonbonded");
  my @libraries = ();
  my @tests = ("psimrcc-sp1", "mcscf-rhf1", "mcscf-rohf1", "mcscf-twocon1", "scf+d-opt1");
  
  foreach my $bin (@binaries) {
    system("rm -rf $dirname/src/bin/$bin");
    system("sed s/$bin// < $dirname/src/bin/Makefile.in > tmp.in");
    system("mv tmp.in $dirname/src/bin/Makefile.in");
    # assuming that each entry in AC_CONFIG_FILES is on its own line
    system("cat $dirname/configure.ac | grep -v 'src/bin/$bin/Makefile' > tmp.ac");
    system("mv tmp.ac $dirname/configure.ac");
  }

  foreach my $lib (@libraries) {
    system("rm -rf $dirname/src/lib/$lib");
    system("sed s/$lib// < $dirname/src/lib/Makefile.in > tmp.in");
    system("mv tmp.in $dirname/src/lib/Makefile.in");
    # assuming that each entry in AC_CONFIG_FILES is on its own line
    system("cat $dirname/configure.ac | grep -v 'src/lib/$lib/Makefile' > tmp.ac");
    system("mv tmp.ac $dirname/configure.ac");
  }

  foreach my $test (@tests) {
    system("rm -rf $dirname/tests/$test");
    system("sed s/$test// < $dirname/tests/Makefile.in > tmp.in");
    system("mv tmp.in $dirname/tests/Makefile.in");
    # assuming that each entry in AC_CONFIG_FILES is on its own line
    system("cat $dirname/configure.ac | grep -v 'tests/$test/Makefile' > tmp.ac");
    system("mv tmp.ac $dirname/configure.ac");
  }

}
