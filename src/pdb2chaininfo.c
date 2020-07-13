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
 * This code listing the user specified identifier-ids from the protein.
 * How to run:
 * ./pdb2chaininfo.exe <PDB-ID> <identifier-id>
 * say,
 * ./pdb2chaininfo.exe 1NHZ.pdb 486
*/
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
int main(int argc, char *argv[])
{
FILE *fr;

if (argc<3) {
  printf("                 <Input file> <Input argument>\n");	
  printf("./pdb2chaininfo.exe <PDB-ID> <identifier-id>\n");
  exit(1);
}

char s[100], res[10];
char hetid[4];
fr=fopen(argv[1],"r");
strcpy(hetid, argv[2]);
//printf("%s\n", hetid);

while(fgets(s,100,fr)!=NULL)
	{
		if(strncmp(s,"HETATM",6)==0 || strncmp(s,"ATOM",4)==0)
		{
		if(strcasestr(s,"ala")==NULL && strcasestr(s,"arg")==NULL && 
		strcasestr(s,"asn")==NULL && strcasestr(s,"asp")==NULL &&
		strcasestr(s,"ash")==NULL && strcasestr(s,"cys")==NULL && 
		strcasestr(s,"cyx")==NULL && strcasestr(s,"glu")==NULL &&
		strcasestr(s,"gln")==NULL && strcasestr(s,"glh")==NULL &&
		strcasestr(s,"gly")==NULL && strcasestr(s,"hie")==NULL &&
		strcasestr(s,"hid")==NULL && strcasestr(s,"hip")==NULL &&
		strcasestr(s,"ile")==NULL && strcasestr(s,"leu")==NULL &&
		strcasestr(s,"lys")==NULL && strcasestr(s,"lyn")==NULL &&
		strcasestr(s,"met")==NULL && strcasestr(s,"phe")==NULL &&
		strcasestr(s,"pro")==NULL && strcasestr(s,"ser")==NULL &&
		strcasestr(s,"thr")==NULL && strcasestr(s,"trp")==NULL &&
		strcasestr(s,"tyr")==NULL && strcasestr(s,"val")==NULL &&
		strcasestr(s,"cym")==NULL && strcasestr(s,"his")==NULL &&
		strcasestr(s,"hoh")==NULL && strcasestr(s,"wat")==NULL)
			{
			res[0]=s[17];
			res[1]=s[18];
			res[2]=s[19];
			res[3]=s[20];
			res[4]=s[21];
			res[5]=s[22];
			res[6]=s[23];
			res[7]=s[24];
			res[8]=s[25];
			res[9]=s[26];
			res[10]='\0';

			if(strstr(res, hetid))
			printf("%s\n",res);
			}
		}
	}
fclose(fr);
}
