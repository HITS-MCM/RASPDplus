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
 * Purpose: Check whether molar refractivity or, wiener or, molecular weight value is zero or, non zero. This program will then select only non zero cases and the corresponding line numbers. 
 */

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
int main(int argc, char *argv[])
{
FILE *fp, *fw1, *fw2;

if (argc<2) {
  printf("./param_chek_nonzero.exe <Input file> <Output file-1> <Output file-2>\n");
  exit(1);
}

char st[500];
int line=1, count=0;
char list[50];
float maxd=0;	
float logp=0, mr=0, wiener=0, molW=0;
int hbd=0, hba=0;
fp=fopen(argv[1],"r");
fw1=fopen(argv[2],"w");
fw2=fopen(argv[3],"w");

if (fw1==NULL || fw2==NULL)            
        {
        printf("param_chek_nonzero.exe: Cannot write output. No output file name is given\n");
        exit(1);
        }

while(fgets(st,500,fp)!=NULL)
{
sscanf(st,"%s %d %d %f %f %f %f %f", list, &hbd, &hba, &logp, &mr, &wiener, &molW, &maxd);

	if(mr>0 && wiener>0 && molW>0)
		{	
			printf("%s\n", list);
			fprintf(fw1,"%d\t %d\t %0.3f\t %0.3f\t %0.3f\t %0.3f\n",hbd, hba, logp, mr, wiener, molW);
			fprintf(fw2, "%0.3f\t%0.3f\n", maxd, maxd);
		}	
}
fclose(fp);
fclose(fw1);
fclose(fw2);
}
