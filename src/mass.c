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
 * Purpose: Calculate moleular weight of small drug-like molecule
 * Input is a pdb file without chain information of a ligand
 */

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include <math.h>
int main(int argc, char *argv[])
{
FILE *fp;
if (argc<2) {
  printf("./mass.exe <Input pdb file>\n");
  exit(1);
}

float sum=0;
char str[100];
float x2,y2,z2;
float m;

char atom[10],t[6],d[6];
int atomno,r;
fp=fopen(argv[1],"r");
int ln=1;
while(fgets(str,100,fp)!=NULL)
	{
	if(strncmp(str,"ATOM",4)==0 || strncmp(str,"HETATM",6)==0)
		{
		sscanf(str,"%s%d%s%s%d%f%f%f",atom,&atomno,t,d,&r,&x2,&y2,&z2);
		if(t[0]!='1' && t[0]!='2' && t[0]!='3' && t[0]!='4' && t[0]!='5' && t[0]!='6' && t[0]!='7' && t[0]!='8' && t[0]!='9' && t[0]!='0')
			{
			if(t[0]=='C'&& t[1]!='l') 
			m=12;
			else if (t[0]=='N' && t[1]!='a')
			m=14;
			else if (t[0]=='O')
			m=16;
			else if (t[0]=='S')
			m=32;
			else if (t[0]=='H')
			m=1;
			if(t[0]=='C'&&t[1]=='l')
			m=35.5;
			if(t[0]=='C'&&t[1]=='L')
			m=35.5;
			else if (t[0]=='P')
			m=31;
			else if (t[0]=='F')
			m=19;
			else if (t[0]=='B')
			m=79;
			else if (t[0]=='I')
			m=126.9;
			else if (t[0]=='N' && t[1]=='a')
			m=23;
			else if (t[0]=='K')
			m=39;
			else if (t[0]=='B')
                        m=10.81;
			else if(t[0]=='S'&&t[1]=='e')
                        m=78.96;
			else if(t[0]=='S'&&t[1]=='E')
                        m=78.96;
			else if(t[0]=='T'&&t[1]=='I')
                        m=47.867;
                        else if(t[0]=='T'&&t[1]=='i')
                        m=47.867;
			else if(t[0]=='S'&&t[1]=='I')
                        m=28.0855;
                        else if(t[0]=='S'&&t[1]=='i')
                        m=28.0855;
			}

		if(t[0]=='1' || t[0]=='2' || t[0]=='3' || t[0]=='4' || t[0]=='5' || t[0]=='6' || t[0]=='7' || t[0]=='8' || t[0]=='9' || t[0]=='0')
			{
			if(t[1]=='C'&& t[2]!='l')
			m=12;
			else if (t[1]=='N' && t[2]!='a')
			m=14;
			else if (t[1]=='O')
			m=16;
			else if (t[1]=='S')
			m=32;
			else if (t[1]=='H')
			m=1;
			if(t[1]=='C'&&t[2]=='l')
			m=35.5;
			if(t[1]=='C'&&t[2]=='L')
			m=35.5;
			else if (t[1]=='P')
			m=31;
			else if (t[1]=='F')
			m=19;
			else if (t[1]=='B')
			m=79;
			else if (t[1]=='I')
			m=126.9;
			else if (t[1]=='N' && t[2]=='a')
			m=23;
			else if (t[1]=='K')
			m=39;
			else if (t[1]=='B')
                        m=10.81;
			if(t[1]=='S'&&t[2]=='e')
                        m=78.96;
                        if(t[1]=='S'&&t[2]=='E')
                        m=78.96;
			else if(t[1]=='T'&&t[2]=='I')
                        m=47.867;
                        else if(t[1]=='T'&&t[2]=='i')
                        m=47.867;
			else if(t[1]=='S'&&t[2]=='I')
                        m=28.0855;
                        else if(t[1]=='S'&&t[2]=='i')
                        m=28.0855;
			}
		sum=sum+m;
		ln++;
		}
	}
printf("%0.1f\n",sum);
fclose(fp);
}
