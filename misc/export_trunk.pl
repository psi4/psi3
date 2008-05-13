#!/usr/bin/perl

#
# Prepare PSI code for shipment to non-developer parties
# Edward Valeev (evaleev at vt.edu) 2008
#

my $svnurl = "https://sirius.chem.vt.edu/svn/psi3/trunk";
my $debug = 0;

my $nargs = $#ARGV;
usage() if ( $nargs < 0 || $nargs > 1 );

my $make_patch = ( $nargs == 1 );
my $rev1       = shift;
my $rev2       = shift if ($make_patch);

# if asked for HEAD, determine the actual revision number from the repository
if ( $rev1 eq "HEAD" ) {
	$_ = `svn info $svnurl | grep 'Revision:'`;
	chomp;
	my @tokens = split;
	$rev1 = $tokens[1];
}

die "REV1 <= REV2" if ( $rev1 <= $rev2 );

my $dirname1 = checkout($rev1);
if ($make_patch) {
    my $dirname2 = checkout($rev2);
    make_patch($dirname1,$dirname2);
}
else {
	make_tgz($dirname1);
}

exit;

0;

sub usage() {

	printf STDOUT "Usage: export_trunk.pl REV1 [REV2]\n";
	printf STDOUT
"where: REV1            -- the target revision to export (numeric revision specifier or HEAD)\n";
	printf STDOUT
"       REV2 (optional) -- if specified, create a patch file that brings REV2 to REV1 (REV1 must be greater than REV2)\n";

}

sub checkout() {
	my $rev = shift;
	
	my $briefdirname = "psi3-r$rev";
	system("rm -rf $briefdirname");
	system("svn export -r $rev $svnurl $briefdirname");

	# extract psi version and buildid from checked-out configure.ac
	my $version;
	my $buildid;
	{
		$_ = `grep psi_version $briefdirname/configure.ac | grep define`;
		chomp;
		s/\[/ /g;
		s/\]/ /g;
		my @tokens = split;
		$version = $tokens[3];
		$_       = `grep psi_buildid $briefdirname/configure.ac | grep define`;
		chomp;
		s/\[/ /g;
		s/\]/ /g;
		@tokens  = split;
		$buildid = $tokens[3];
	}

    print "Working on: psi3 version $version buildid $buildid revision $rev\n";

	# form the full name and rename the directory
	my $dirname = "psi-$version-$buildid-r$rev";
	system("mv $briefdirname $dirname");

	# sanitize
	sanitize($dirname);

	# prepare for build
	system("(cd $dirname; autoconf; rm -rf autom4te.cache)");

	return $dirname;
}

sub make_tgz() {
	my $dirname = shift;
	# make .tgz
    system("rm -rf $dirname.tar.gz");
    system("tar cvzf $dirname.tar.gz $dirname");
    system("rm -rf $dirname");
    print "Prepared $dirname.tar.gz\n";
}

sub make_patch() {
    my $dirname1 = shift;
    my $dirname2 = shift;
    
    my ($version1,$buildid1,$rev1) = parse_dirname($dirname1);
    my ($version2,$buildid2,$rev2) = parse_dirname($dirname2);
    my $patchname = "psi-patch-$rev2-to-$rev1";
    system("rm -rf $patchname");
    system("diff -Naur $dirname2 $dirname1 > $patchname");
    system("rm -rf $dirname1 $dirname2");
    print "Prepared $patchname\n";
}

sub parse_dirname() {
	my $dirname = shift;
	my ($psi, $version, $buildid, $rrev) = split /-/, $dirname;
	
	return ($version, $buildid, $rrev);
}

#
# Rip out pieces not ready for distribution
#
sub sanitize() {

    my $dirname = shift;
    
	system("rm -rf $dirname/misc");

	my @binaries  = ( "psimrcc", "mcscf", "nonbonded" );
	my @libraries = ();
	my @tests     = (
		"psimrcc-sp1", "mcscf-rhf1", "mcscf-rohf1", "mcscf-twocon1",
		"scf+d-opt1"
	);

	foreach my $bin (@binaries) {
		system("rm -rf $dirname/src/bin/$bin");
		system("sed s/$bin// < $dirname/src/bin/Makefile.in > tmp.in");
		system("mv tmp.in $dirname/src/bin/Makefile.in");

		# assuming that each entry in AC_CONFIG_FILES is on its own line
		system(
"cat $dirname/configure.ac | grep -v 'src/bin/$bin/Makefile' > tmp.ac"
		);
		system("mv tmp.ac $dirname/configure.ac");
	}

	foreach my $lib (@libraries) {
		system("rm -rf $dirname/src/lib/$lib");
		system("sed s/$lib// < $dirname/src/lib/Makefile.in > tmp.in");
		system("mv tmp.in $dirname/src/lib/Makefile.in");

		# assuming that each entry in AC_CONFIG_FILES is on its own line
		system(
"cat $dirname/configure.ac | grep -v 'src/lib/$lib/Makefile' > tmp.ac"
		);
		system("mv tmp.ac $dirname/configure.ac");
	}

	foreach my $test (@tests) {
		system("rm -rf $dirname/tests/$test");
		system("sed s/$test// < $dirname/tests/Makefile.in > tmp.in");
		system("mv tmp.in $dirname/tests/Makefile.in");

		# assuming that each entry in AC_CONFIG_FILES is on its own line
		system(
"cat $dirname/configure.ac | grep -v 'tests/$test/Makefile' > tmp.ac"
		);
		system("mv tmp.ac $dirname/configure.ac");
	}

}
