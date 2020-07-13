#!/usr/bin/perl

#>  Copyright (c) 2013 2014 2015 2016 2017 2018 2019 2020
#>  Heidelberg University and Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org)
#>  Schloss-Wolfsbrunnenweg 35
#>  69118 Heidelberg, Germany
#>
#>  Supercomputing Facility for Bioinformatics and Computational Biology, IIT Delhi (http://www.scfbio-iitd.res.in/)
#>  Indian Institute of Technology
#>  Hauz Khas, New Delhi - 110016, India
#>
#>  Please send your contact address to get information on updates and
#>  new features to "mcmsoft@h-its.org". Questions will be
#>  answered as soon as possible.
#>  References:
#>  A rapid identification of hit molecules for target proteins via physico-chemical descriptors.
#>  (2013) Phys. Chem. Chem. Phys., 15, 9107-9116.
#>  Authors: Goutam Mukherjee and B. Jayaram
#>  Version 1.0 (April 2013)

#>  RASPD+: Fast protein-ligand binding free energy prediction using simplified physicochemical features
#>  Authors: Stefan Holderbach, Lukas Adam, B. Jayaram, Rebecca C. Wade, Goutam Mukherjee
#>  Version 1.0 (June 2020)
#>  ChemRxiv preprint (https://doi.org/10.26434/chemrxiv.12636704.v1), 2020

my $atomline=0;

@STORE;

open(F1,$ARGV[0])|| die("Can't open file:".$ARGV[0]);

while(<F1>)
			{
				  $atomline++  if(/^ATOM/); 
    				  next if(!/^0/);
   				  @lcontents = split(/\s+/,$_);
				  ($val) = @lcontents;
				  
				  for( $i = 0; $i < length $val; $i+=3) {
					
				   					 $sub = substr ( $val , $i, 3);
				   					 $sub =~s/^0//g;
				   					 $sub =~s/^0//g;
				   
				   					 $temp1 = $sub if($i == 0);
				   					 $temp2 = $sub if($i == 3);
				   
				   if($i == 6) {
							#print "$temp1 $temp2 $sub\n";
							$STORE[$temp1] += $sub;
							$STORE[$temp2] += $sub;
							#print "$STORE[$temp1] $STORE[$temp2]\n";
						}
									 }

			}
close(F1);
				  #print "BOND ORDER\n";
				  for( $i = 1; $i <= $atomline; $i++) {
								#print "$i  =>  $STORE[$i] \n";
								print "$STORE[$i] \n";
								}
