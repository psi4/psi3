#!/usr/bin/perl

use Cwd;

$distdir = cwd();
$version = "PSI 3.2";

system ("clear");
print "***********************************\n";
print "*   $version Installation Script   *\n";
print "***********************************\n\n";

print "STEP 1:\n";
print "Running autoconf ...\n\n";
system ("autoconf");

print "STEP 2:\n";
print "Where would you like to install the $version binaries?\n";
print "We recommend the default path /usr/local/psi. If the\n";
print "default path is correct just hit enter, if not please\n";
print "type the full path:\n";
chomp($prefix = <STDIN>);
print "\n";

$arch = `./bin/host.sh`;
chomp($arch);

print "STEP 3:\n";
print "Where would you like the $version compilation directory?\n";
print "If you need executables for several architectures, it may\n";
print "be helpful to name the compilation directory by the \n";
print "architecture type, which is $arch for the present machine:\n";
chomp($compdir = <STDIN>);
print "\n";

print "STEP 4:\n";
print "Checking for BLAS and LAPACK ...\n";
print "\n";
$dir = "/usr/lib";
@lapack = `ls $dir/*lapack*`;

$found = 0;
foreach $file (@lapack) {
  chomp($file);
  if ($file eq "$dir/liblapack.a") {
    $found = 1;
  }
}

if ($found == 1) {
  print "Found the LAPACK library in $dir\n\n";
  $lapack = "-llapack";
}
else {
  print "Did not find the LAPACK library\n";
  print "Please enter the location of the LAPACK library:\n";
  chomp($lapack = <STDIN>);
  print "\n";
}

@blas = `ls $dir/*blas*`;
$found = 0;
foreach $file (@blas) {
  chomp ($file);
  if ($file eq "$dir/libblas.a") {
    $found = 1;
  }
}

if ($found == 1) {
  print "Found the BLAS library in $dir\n\n";
  $blas = "-lblas";
}
else {
  print "Did not find the BLAS library\n";
  print "Please enter the location of the BLAS library:\n";
  chomp($blas = <STDIN>);
  print "\n";
}

print "STEP 5:\n";
print "Making the $version compilation directory ...\n\n";
system ("mkdir $compdir");

print "Moving to the $version compilation directory ...\n\n";
chdir ("$compdir");

print "STEP 6:\n";
print "Running the configure script ... \n\n";
system ("$distdir/configure --prefix=$prefix -with-blas=$blas -with-lapack=$lapack");

print "STEP 7:\n";
print "Compiling $version ...\n";
system ("make");

system ("clear");

print "STEP 8:\n";
TEST: {
  print "Would you like to run the $version test cases? (Y|N)\n";
  chomp ($test = <STDIN>);
  redo TEST unless $test eq "Y" || "N";
}
if ($test eq "Y") {
  system ("make tests");
}
elsif ($test eq "N") {
  break;
}

print "STEP 9:\n";
print "Installing $version ...\n";
system ("make install");

print "\n\n";
print "**************************************\n";
print "*   $version Installation Complete    *\n";
print "**************************************\n";
print "\n\n";

