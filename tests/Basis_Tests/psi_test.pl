#!/usr/bin/perl 
########################################################################
#
#	This script takes a standard PSI input and performs
#	atomic calculations on the H-Ar specifiying the proper 
#	multiplicity, occupation, reference ... etc for H-Ar.
#
#	Options - 0 - will perform test on ~4o basis sets given below
#		    ( the starting input.dat should be identical to the
#		    file named standart.input )
#		- 1 - will perform test on basis set given in input.dat
#		    ( input file must start with hydrogen )
#
#	Note: the input file should always start wit
#
#	It writes the atomic energies into a file named 
#	<reference>.<basis set>.Energies 
#
#
#
#	Dec. 2002
########################################################################

# HASH OF ELEMENT NAMES

%elemname = (
  h  => HYDROGEN,
  he => HELIUM,
  li => LITHIUM,
  be => BERYLLIUM,
  b  => BORON,
  c  => CARBON,
  n  => NITROGEN,
  o  => OXYGEN,
  f  => FLUORINE,
  ne => NEON,
  na => SODIUM,
  mg => MAGNESIUM,
  al => ALUMINUM,
  si => SILICON,
  p  => PHOSPHORUS,
  s  => SULFUR,
  cl => CHLORINE,
  ar  => ARGON,
  h => HYDROGEN);

# HASH OF ELEMENT MULTIPLICITIES

%elemult = (
  h  => 2,
  he => 1,
  li => 2,
  be => 1,
  b  => 2,
  c  => 3,
  n  => 4,
  o  => 3,
  f  => 2,
  ne => 1,
  na => 2,
  mg => 1,
  al => 2,
  si => 3,
  p  => 4,
  s  => 3,
  cl => 2,
  ar => 1,
  h => 2);

# ARRAY OF ELEMENT SYMBOLS 

@elemArray = qw(h he li be b c n o f ne na mg al si p s cl ar h);

# ARRAY OF BASIS SETS

@Basis= ("DZ", "DZP", "DZ-DIF", "DZP-DIF",
        "WACHTERS", "WACHTERS-F",
        "STO-3G", "3-21G", "6-31G","6-311G",
        "TZ2P", "TZ2PF", "TZ-DIF", "TZ2P-DIF", "TZ2PF-DIF",
        "CC-PVDZ", "CC-PVTZ", "CC-PVQZ", "CC-PV5Z", "CC-PV6Z",
        "CC-PCVDZ", "CC-PCVTZ", "CC-PCVQZ", "CC-PCV5Z", "CC-PCV6Z",
        "AUG-CC-PVDZ", "AUG-CC-PVTZ", "AUG-CC-PVQZ", "AUG-CC-PV5Z", "AUG-CC-PV6Z",
        "AUG-CC-PCVDZ", "AUG-CC-PCVTZ", "AUG-CC-PCVQZ", "AUG-CC-PCV5Z", "AUG-CC-PCV6Z",
        "PV7Z", "AUG-PV7Z", "AUG-CC-PV7Z");

# CHOOSE BETWEEN A SINGLE BASIS SET AND ALL BASIS SETS

print("Enter 0 - to do a complete test case\n");
print("      1 - to test a specific basis set\n");
$choice = <STDIN>;

# LOOP OVER BASIS SETS AND ELEMENTS

for($bas=0; $bas<40; $bas++){
  for($i=0;$i<19;$i++)
  {
     $energy[$i]=0;
     Read_input();
     system("input");
     system("psi"); 
     $energy[$i]=0;
     energy($energy[$i]); 
     print "Energy = $energy[$i]\n"; 
     system("mv output.dat $elemArray[$i].out"); 
     system("psiclean"); 
  }
  print_data();
  print "The wavefunction is $wfn \n";
  print "The basis set is $basis\n puream = $puream\n";
  rename("data", "$wfn.$basis.Energies");

# TRUNCATION
  if(($choice==1) || ($Basis[$bas]=/AUG-CC-PV7Z/))
	{$bas=55;}
}
###############################################################

# READS AND MODIFIES INPUT FILE

sub Read_input 
{
  #OPEN INPUT.DAT FOR PSI3.0
  open(IN, "+<input.dat") || die "cannot open input.dat: $!";

# $ename = $elemArray[$i-1]; 
  $count=0;
  @input = <IN>;
  seek(IN,0,0);
  while(<IN>)
  {
    $count++;
    
   # READ TYPE OF WAVEFUNCTION
   if($input[$count]=~/wfn/)
    {
	$tmp = $input[$count];
	@tmp2 = split(/ +/, $tmp);
	chomp($tmp2[3]);
	$wfn = $tmp2[3];
    }	

    # MODIFY REFERENCE TYPE
    if($input[$count]=~/reference/)
    {
	if($elemult{$elemArray[$i]} == 1)
	{
        $tmp = $input[$count];
        @tmp2 = split(/ +/, $tmp);
        chomp($tmp2[3]);
	$input[$count] =~ s/$tmp2[3]/rhf/;
       $reference = "rhf";
	}
	elsif($elemult{$elemArray[$i]}!=1)
       {
        $tmp = $input[$count];
        @tmp2 = split(/ +/, $tmp);
        chomp($tmp2[3]);
	$input[$count] =~ s/$tmp2[3]/rohf/;
       $reference = "uhf";
        }

    }

   # CHANGE MULTIPLICITY
   if($input[$count]=~/multp/)
    {
    $input[$count] =~ s/[0-9]/$elemult{$elemArray[$i]}/;
    }

    # CHANGE ATOM
    if($input[$count]=~/zmat/)
    {
    $input[$count+1] =~ s/$elemArray[$i-1]/$elemArray[$i]/;
    }

   # SPECIFY PROPER GROUND STATE OCCUPATION
   # first the DOCC   
     if($input[$count]=~/docc/)
	{
        $tmp = $input[$count];
        @tmp2 = split(/ = /, $tmp);
        chomp($tmp2[1]);
        if(($i==0) || ($i==18)) 
        	{ $input[$count] =~ s/$tmp2[1]/0 0 0 0 0 0 0 0/;}
        if(($i==1) || ($i ==2)) 
                { $input[$count] =~ s/$tmp2[1]/1 0 0 0 0 0 0 0/;}
        if(($i>2) && ($i<7)) 
                { $input[$count] =~ s/$tmp2[1]/2 0 0 0 0 0 0 0/;}
        if($i==7)
                { $input[$count] =~ s/$tmp2[1]/2 0 0 0 0 1 0 0/;}
        if($i==8)    
                { $input[$count] =~ s/$tmp2[1]/2 0 0 0 0 1 1 0/;}
	if(($i==9) || ($i==10))
		{ $input[$count] =~ s/$tmp2[1]/2 0 0 0 0 1 1 1/;}
        if(($i>10) && ($i<15))
                { $input[$count] =~ s/$tmp2[1]/3 0 0 0 0 1 1 1/;}
        if($i==15)
                { $input[$count] =~ s/$tmp2[1]/3 0 0 0 0 2 1 1/;}
        if($i==16)
                { $input[$count] =~ s/$tmp2[1]/3 0 0 0 0 2 2 1/;}
        if($i==17)
                { $input[$count] =~ s/$tmp2[1]/3 0 0 0 0 2 2 2/;}
	}

	# Then the SOCC
       if($input[$count]=~/socc/)
        {
        $tmp = $input[$count];
        @tmp2 = split(/ = /, $tmp);
        chomp($tmp2[1]);
        if(($i==0) || ($i==18))             
                { $input[$count] =~ s/$tmp2[1]/1 0 0 0 0 0 0 0/;}
        if((($i==1)||($i==9)) || ($i==17))    
                { $input[$count] =~ s/$tmp2[1]/0 0 0 0 0 0 0 0/;}
        if(($i==2) || ($i==10))    
                { $input[$count] =~ s/$tmp2[1]/1 0 0 0 0 0 0 0/;}
        if(($i==3) || ($i==11))  
                { $input[$count] =~ s/$tmp2[1]/0 0 0 0 0 0 0 0/;}
        if(($i==4) || ($i==12))  
                { $input[$count] =~ s/$tmp2[1]/0 0 0 0 0 1 0 0/;}
        if(($i==5) || ($i==13))  
                { $input[$count] =~ s/$tmp2[1]/0 0 0 0 0 1 1 0/;}
        if(($i==6) || ($i==14))
                { $input[$count] =~ s/$tmp2[1]/0 0 0 0 0 1 1 1/;}
        if(($i==7) || ($i==15)) 
                { $input[$count] =~ s/$tmp2[1]/0 0 0 0 0 0 1 1/;}
        if(($i==8) || ($i==16))
                { $input[$count] =~ s/$tmp2[1]/0 0 0 0 0 0 0 1/;}
        }                      
	
		 		
    # READ ANR MODIFY BASIS SET NAME
      if($input[$count]=~/basis/)
        {
        $tmp3 = $input[$count];
        @tmp4 = split(/"/, $tmp3);
	if($choice==1)
        	{$basis = $tmp4[1];}
	else{
        	$basis = $Basis[$bas];
        	$input[$count] =~ s/$tmp4[1]/$basis/;
		}	
        }
    
	# CHECK IF puream = true
        if($input[$count]=~/puream/)
        {
        	if(($basis=~/cc\b/i) || ($basis=~/6-311/))
	  	{
	        $tmp = $input[$count];
       		@tmp2 = split(/ +/, $tmp);
        	chomp($tmp2[3]);
        	$input[$count] =~ s/$tmp2[3]/true/;
		$puream  = "true";
	  	}
		else
		{
                $tmp = $input[$count];
                @tmp2 = split(/ +/, $tmp);
                chomp($tmp2[3]);
                $input[$count] =~ s/$tmp2[3]/false/;
                $puream  = "false";
		}
	}
  } 
  seek(IN,0,0);
  print IN @input;
  truncate(IN,tell(IN));
  close (IN);
}
###############################################################

# READ ENERGY FROM OUTPUT FILE

sub energy
{
open(OUT, "output.dat") || die "cannot open output.dat: $!";

while (<OUT>)
  {
# IF USING DETCI
#    if($wfn = ~/fci/)
#	{
#		if(/ROOT 1/)
#    		{
#      		$line = $_;
#      		@data = split(/ +/,$line);
#      		chomp($data[4]);
#      		$energy[$i] = $data[4];
#    		}
#	}
#    if($wfn = ~/scf/)
#	{
#  IF WFN IS SCF
		if(/total energy/)
		{
       		$line = $_;	
		@junk = split(/ +/, $line);
		chomp($junk[4]);
        	$_[0] = $junk[4];
		}
# 	}
  }
close(OUT);
}
#################################################################

# PRINT DATA TO FILE NAMED <reference>.<basis set>.Energies
sub print_data
{
open(DAT, ">data") || die "cannot open data: $!";

 print DAT "\t\t	$basis\n\n";
 print DAT "  ATOM\t  	   MULTIPLICITY\t       ENERGY\n";
 print DAT "___________\t  _______________\t  _________________\n";

format DAT = 
@<<<<<<<<<<<<<<<<<< @<<<<<<<<<<<<<<<<<< @<<<<<<<<<<<<<<<<<<< 
$elemname{$elemArray[$i]},        $elemult{$elemArray[$i]},        $energy[$i],
.

for($i=0; $i<18; $i++)
	{
	write(DAT);
	}
close(DAT);
}

###################	END	################################
