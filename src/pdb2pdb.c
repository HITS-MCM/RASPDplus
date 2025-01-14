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
 * Purpose: Format a pdf file of ligand (non-standard) molecule.
 * How to run:
 * ./pdb2pdb.exe <ligand name with extension (.pdb)> 
 * say,
 * ./pdb2pdb.exe ligand.pdb
*/
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
int main(int argc, char *argv[])
{
	FILE *fr;
	if (argc<2) {
        printf("./pdb2pdb.exe <Input file>\n");
        exit(1);
        }
	char s[100];
	int resno=1, atomno=1;

	char atom[6]="ATOM", atmsmb[4];
	char resid[4]="DRG";
	char xcod[8],ycod[8],zcod[8];
	float xcod1,ycod1,zcod1;

	fr=fopen(argv[1],"r");

	while(fgets(s,100,fr)!=NULL)
		{
		if(strncmp(s,"ATOM",4)==0||strncmp(s,"HETATM",6)==0)
			{
			atmsmb[0]=s[12];
			atmsmb[1]=s[13];
			atmsmb[2]=s[14];
			atmsmb[3]='\0';

			xcod[0]=s[30];
			xcod[1]=s[31];
			xcod[2]=s[32];
			xcod[3]=s[33];
			xcod[4]=s[34];
			xcod[5]=s[35];
			xcod[6]=s[36];
			xcod[7]=s[37];
			xcod[8]='\0';
			xcod1=atof(xcod);

			ycod[0]=s[38];
			ycod[1]=s[39];
			ycod[2]=s[40];
			ycod[3]=s[41];
			ycod[4]=s[42];
			ycod[5]=s[43];
			ycod[6]=s[44];
			ycod[7]=s[45];
			ycod[8]='\0';
			ycod1=atof(ycod);

			zcod[0]=s[46];
			zcod[1]=s[47];
			zcod[2]=s[48];
			zcod[3]=s[49];
			zcod[4]=s[50];
			zcod[5]=s[51];
			zcod[6]=s[52];
			zcod[7]=s[53];
			zcod[8]='\0';
			zcod1=atof(zcod);
				
			strncpy(resid, "DRG", 3);
			printf("%-6s%5d %-4s %3s   %3d    %8.3f%8.3f%8.3f\n",atom,atomno,atmsmb,resid,resno,xcod1,ycod1,zcod1); 
			atomno++;
			}
		}

fclose(fr);
}

