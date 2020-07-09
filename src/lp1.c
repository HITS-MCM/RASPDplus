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
#>  Manuscript in preparation, 2020

 * Written by Goutam Mukherjee
 *Purpose: Print atom symbol of atoms present in a ligand (non-standard) molecule.
 * How to run:
 * ./lp1.exe <ligand name with extension (.pdb)> <output file>
 * say,
 * ./pdb2pdb.exe ligand.pdb atomsymbol
*/

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
//#include <iostream>
//#include <cstdio>
//#include <cstring>
//#include <vector>
//using namespace std;

int main(int argc, char *argv[])
{
FILE *fp,*fw;
struct format
{
  char atom[6],atmsmb[3],resid[4];
  int atomno,resno,k;
  float x,y,z;
};
struct format f;

char s[500];

if (argc<2) {
  printf("./lp1.exe <ligand name with extension (.pdb)> <output file>\n");
  exit(1);
}

fp=fopen(argv[1],"r");
fw=fopen(argv[2],"w");


if (fw==NULL) 
	{
  	printf("lp1.exe: Cannot write output. No output file name is given\n");
  	exit(1);
	}


while(fgets(s,500,fp)!=NULL)
	{
	if(strncmp("ATOM",s,4)==0||strncmp("HETATM",s,6)==0)
		{
		sscanf(s,"%s%d%s%s%d%f%f%f",f.atom, &f.atomno, f.atmsmb, f.resid, &f.resno, &f.x, &f.y, &f.z);

		if((f.atmsmb[0]=='1')||(f.atmsmb[0]=='2')||(f.atmsmb[0]=='3')||(f.atmsmb[0]=='4')||(f.atmsmb[0]=='5')||(f.atmsmb[0]=='6')||(f.atmsmb[0]=='7')||(f.atmsmb[0]=='8')||(f.atmsmb[0]=='9')||(f.atmsmb[0]=='0')||strncmp(f.atmsmb,"'",1)==0||strncmp(f.atmsmb,"*",1)==0)
			{
			if(strstr(f.atmsmb,"Cl")==NULL)
			fprintf(fw,"%c",f.atmsmb[1]);
			if(strstr(f.atmsmb,"Cl"))
			fprintf(fw,"%c",f.atmsmb[2]);
			}
		if((f.atmsmb[0]!='1')&&(f.atmsmb[0]!='2')&&(f.atmsmb[0]!='3')&&(f.atmsmb[0]!='4')&&(f.atmsmb[0]!='5')&&(f.atmsmb[0]!='6')&&(f.atmsmb[0]!='7')&&(f.atmsmb[0]!='8')&&(f.atmsmb[0]!='9')&&(f.atmsmb[0]!='0')&&(f.atmsmb[0]!='0')&&strncmp(f.atmsmb,"'",1)!=0&&strncmp(f.atmsmb,"*",1)!=0)
			{
			if(strstr(f.atmsmb,"Cl")==NULL)
			fprintf(fw,"%c",f.atmsmb[0]);
			if(strstr(f.atmsmb,"Cl"))
			fprintf(fw,"%c",f.atmsmb[1]);
			}
		}
	}
fclose(fp);
}

