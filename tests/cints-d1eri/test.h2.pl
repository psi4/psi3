#!/usr/bin/perl

#
# tests output of cints --deriv1_ints when run on H2 molecule
#

die if ($#ARGV != 0);

my $output_dat = $ARGV[0];

# parse integrals first
open IFILE, "<$output_dat" || die "could not open $output_dat";
my $cints_count = 0;
while (<IFILE>) {
  $cints_count++ if /CINTS/;
  next if $cints_count<3;

  # find the integrals
  goto PARSE_SO_DERIV if /Electron repulsion integrals/;
}

PARSE_SO_DERIV:
my %integrals;
<IFILE>;
while (<IFILE>) {
  goto DONE_PARSING_INTS if ($_ eq "\n");
  my ($der, $i, $j, $k, $l, $val) = split;
  #printf STDOUT "Value = %lf\n",$val;
  my $key = canon_integral_key($der,$i,$j,$k,$l);
  $integrals{$key} = $val;
}

DONE_PARSING_INTS:
close OFILE;

# now parse MO coefficients
open OFILE, "<$output_dat" || die "could not open $output_dat";
while (<OFILE>) {
  goto SKIPPED_HEADER if /molecular orbitals for irrep Ag/;
}
SKIPPED_HEADER:
# skip 3 lines
<OFILE>; <OFILE>; <OFILE>;
my @mo_coefs;
while (<OFILE>) {
  goto PARSED_MO_COEFS if ($_ eq "\n");
  my ($so, $coef, $rest) = split;
  push @mo_coefs, $coef;
}
PARSED_MO_COEFS:

# finally find two-body skeleton contribution to the gradient
my $scf_te_deriv;
while (<OFILE>) {
  next if !/Two-electron contribution to the forces/;
  # skip next 2 lines
  $_ = <OFILE>; $_ = <OFILE>; $_ = <OFILE>;
  # grab the derivative contribution along z-direction
  my (@tokens) = split;
  $scf_te_deriv = $tokens[3];
  goto PARSED_SCF_TE_DERIV;
}
PARSED_SCF_TE_DERIV:
close OFILE;

my $nso = $#mo_coefs + 1;

my $mo_int = 0.0;
foreach $i (0..$nso-1) {
  foreach $j (0..$nso-1) {
    foreach $k (0..$nso-1) {
      foreach $l (0..$nso-1) {
        my $key = canon_integral_key(0,$i,$j,$k,$l);
        my $so_int = $integrals{$key};
        $mo_int += $mo_coefs[$i] *
                   $mo_coefs[$j] *
                   $mo_coefs[$k] *
                   $mo_coefs[$l] *
                   $so_int;
      }
    }
  }
}

printf STDOUT "scf te deriv along z axis         = %15.10lf\n", $scf_te_deriv;
# symmery-adapted perturbation operator need to be transformed back to nuclear perturbation basis
printf STDOUT "same computed using SO derivative = %15.10lf\n", -1.0*$mo_int/sqrt(2.0);


exit 0;

# compute the canonical key for the integral
sub canon_integral_key
{
  my ($der, $i, $j, $k, $l) = @_;
  my $key;
  
  if ($i < $j) {
    my $tmp = $j;
    $j = $i;
    $i = $tmp;
  }
  if ($k < $l) {
    my $tmp = $l;
    $l = $k;
    $k = $tmp;
  }
  my $ij = $j + $i*($i+1)/2;
  my $kl = $l + $k*($k+1)/2;
  if ($ij < $kl) {
    my $tmp = $k;
    $k = $i;
    $i = $tmp;
    $tmp = $l;
    $l = $j;
    $j = $tmp;
  }
  
  $key = "$der $i $j $k $l";
  
  return $key;
}
