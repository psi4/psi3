#!/usr/bin/perl

use Cwd;

$distdir = cwd();
$version = "PSI 3.2";

system ("clear");

print "***********************************\n";
print "*   $version Installation Script   *\n";
print "***********************************\n\n";

print "STEP 1: Run autoconf\n";
print "STEP 2: Choose directory for binaries\n";
print "STEP 3: Choose compilation directory\n";
print "STEP 4: Locate BLAS and LAPACK\n";
print "STEP 5: Move to the compilation directory\n";
print "STEP 6: Run configure with the correct options\n";
print "STEP 7: Compile $version\n";
print "STEP 8: Run the test suite\n";
print "STEP 9: Install $version\n\n";
print "(Press enter to continue)\n";
$tmp = <STDIN>;

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
print "type the full path:\n\n";
chomp($prefix = <STDIN>);
print "\n";

# Default location for prefix
if ($prefix eq "") {
  $prefix = "/usr/local/psi";
}

print "$version binaries will be located in: $prefix\n\n";

$arch = `./bin/host.sh`;
chomp($arch);

print "STEP 3:\n";
print "Where would you like the $version compilation directory?\n";
print "If you need executables for several architectures, it may\n";
print "be helpful to name the compilation directory by the \n";
print "architecture type, which is $arch for the\n";
print "present machine:\n\n";
chomp($compdir = <STDIN>);
print "\n";

# Default location for compdir
if ($compdir eq "") {
  $compdir = "$distdir/$arch";
}

print "$version compilation directory is: $compdir\n\n";
print "(Press enter to continue)\n";
$tmp = <STDIN>;

system ("clear");

print "***********************************\n";
print "*   $version Installation Script   *\n";
print "***********************************\n\n";

# Search /usr/lib and /usr/local/lib for BLAS and LAPACK

print "STEP 4:\n";
print "Checking for BLAS and LAPACK ...\n\n";

$dir1 = "/usr/lib";
$dir2 = "/usr/local/lib";

@blas = `ls $dir1/*blas*`;
$found = 0;
foreach $file (@blas) {
  chomp ($file);
  if ($file eq "$dir1/libblas.a") {
    $found = 1;
  }
}

if ($found != 1) {
  @blas = `ls $dir2/*blas*`;
  foreach $file (@blas) {
    chomp($file);
    if ($file eq "$dir2/libblas.a") {
      $found = 2;
    }
  }
}

if ($found == 2) {
  $blas = "$dir2/libblas.a";
}
elsif ($found == 1) {
  $blas = "$dir1/libblas.a";
}

if ($found == 2 || $found == 1) {
  print "The default BLAS library is $blas.\n";
  print "If this is not the BLAS library you would like to use,\n";
  print "please enter the full path to the desired BLAS library:\n\n";
  chomp($newblas = <STDIN>);
  print "\n";
  
  if ($newblas ne "") {
    $blas = $newblas;
  }
}
elsif ($found == 0) {
  print "The BLAS library is not located in $dir2 or $dir1.\n";
  print "Please enter the full path to the desired BLAS library:\n\n";
  chomp($blas = <STDIN>);
  print "\n";
}

@lapack = `ls $dir1/*lapack*`;
$found = 0;
foreach $file (@lapack) {
  chomp($file);
  if ($file eq "$dir1/liblapack.a") {
    $found = 1;
  }
}

if ($found != 1) {
  @lapack = `ls $dir2/*lapack*`;
  foreach $file (@lapack) {
    chomp($file);
    if ($file eq "$dir2/liblapack.a") {
      $found = 2;
    }
  }
}

if ($found == 2) {
  $lapack = "$dir2/liblapack.a";
}
elsif ($found == 1) {
  $lapack = "$dir1/liblapack.a";
}

if ($found == 2 || $found == 1) {
  print "The default LAPACK library is $lapack.\n"; 
  print "If this is not the LAPACK library you would like to use,\n";
  print "please enter the full path to the desired LAPACK library:\n\n";
  chomp($newlapack = <STDIN>);
  print "\n";

  if ($newlapack ne "") {
    $lapack = $newlapack;
  }
}
elsif ($found == 0) {
  print "The LAPACK library is not located in $dir2 or $dir1.\n";
  print "Please enter the full path to the desired LAPACK library:\n\n";
  chomp($lapack = <STDIN>);
  print "\n";
}

print "(Press enter to continue)\n";
$tmp = <STDIN>;

system ("clear");

print "***********************************\n";
print "*   $version Installation Script   *\n";
print "***********************************\n\n";

print "STEP 5:\n";
print "Making the $version compilation directory ...\n\n";
system ("mkdir $compdir");

print "Moving to the $version compilation directory ...\n\n";
chdir ("$compdir");

print "STEP 6:\n";
print "The configure script will run with the following options:\n\n";
print "--prefix=$prefix\n";
print "-with-blas='$blas'\n";
print "-with-lapack='$lapack'\n\n";
print "If this is correct press enter.  If this is not correct,\n";
print "please type the options you would like to pass to the\n";
print "configure script:\n\n";
chomp($configure=<STDIN>);

if ($configure eq "") {
  print "Running the configure script ... \n\n";
  system ("$distdir/configure --prefix=$prefix -with-blas='$blas' -with-lapack='$lapack'");
}
else {
  system ("$distdir/configure $configure");
}

print "\n(Press enter to continue)\n";
$tmp = <STDIN>;

system ("clear");

print "***********************************\n";
print "*   $version Installation Script   *\n";
print "***********************************\n\n";

print "STEP 7:\n";
print "Compiling $version ...\n";
print "This could take a while.  You might get a cup of coffee\n";
print "or program CCSDTQ while you wait :)\n\n";
print "(Press enter to begin)\n";
$tmp = <STDIN>;
system ("make");

system ("clear");

print "***********************************\n";
print "*   $version Installation Script   *\n";
print "***********************************\n\n";

print "STEP 8:\n";
TEST: {
  print "Would you like to run the $version test cases? (Y|N)\n\n";
  chomp ($test = <STDIN>);
  redo TEST unless $test eq "Y" || "N";
}
if ($test eq "Y") {
  print "\nRunning test cases ...\n\n";
  print "This could also take a while.  You might get another\n";
  print "cup of coffee or program CCSDTQP while you wait :)\n\n";
  unlink <test.out>;
  system ("make tests > test.out");
  open(IN, "<test.out") || die "cannot open test.out: $!";
  while (<IN>) {
    if (/FAILED/) {
      print "One of the test cases failed!\n";
      print "We do not recomment that you install $version!\n";
      print "Please contact: psimaster@sirius.chem.vt.edu\n";
      exit;
    }
  }
  
  print "The test cases were successful!\n\n";
  print "(Press enter to continue)\n";
  $tmp = <STDIN>;
}
elsif ($test eq "N") {
  break;
}

system ("clear");

print "***********************************\n";
print "*   $version Installation Script   *\n";
print "***********************************\n\n";

print "STEP 9:\n";
print "Installing $version ...\n\n";
print "(Press enter to begin)\n";
$tmp = <STDIN>;
system ("make install");

system ("clear");

print "**************************************\n";
print "*   $version Installation Complete    *\n";
print "**************************************\n";
print "\n\n";

