/*
#>  Copyright (c) 2013 2014 2015 2016 2017 2018 2019 2020
#>  Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org)
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

 * Written by Goutam Mukherjee
 * Purpose: Remove hydrogen atoms from a protein.
 * How to run:
 * ./pdb2noh.exe <PDB-ID> 
 * say,
 * ./pdb2noh.exe 1a30.pdb (source: wget http://www.scfbio-iitd.res.in/software/drugdesign/Method1/1a30.pdb)
*/

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
int main(int argc, char *argv[])
	{
	FILE *fr;
	if (argc<2) {
  	printf("./pdb2noh.exe <Input file>\n");
  	exit(1);
	}
	char s[100];

		struct format
		{	
  		char atom[6], atmsmb[6];
  		int atomno=0;
		};
		struct format f;
fr=fopen(argv[1],"r");

		while(fgets(s, 100, fr)!=NULL)
		{
			if(strncmp(s,"ATOM",4)==0 || strncmp(s,"HETATM",6)==0)
			{
			sscanf(s,"%s %d %s",f.atom, &f.atomno, f.atmsmb);
			if(strncmp(f.atmsmb,"H",1)!=0 && strncmp(f.atmsmb,"1H",2)!=0 && strncmp(f.atmsmb,"2H",2)!=0 && strncmp(f.atmsmb,"3H",2)!=0 && strncmp(f.atmsmb,"4H",2)!=0)
			printf("%s",s);
			}
		}
fclose(fr);
	}

