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
 * Purpose: This program will concatenate hydrogen bond donor, acceptor, logP, molar refractivity and wiener index values in a single file. 
 * How to run:
 * ./lipi.exe <hydrogen bond donor> <hydrogen bond acceptor> <lopP value> <molar refractivity value> <wiener index>
*/

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include <math.h>
int main(int argc, char *argv[])
	{
	FILE *f1,*f2,*f3,*f4,*f5;

	if (argc<6)
        {
        printf("             <Input file-1>           <Input file-2>     <Input file-3>     <Input file-4>        <Input file-5>\n");
        printf("./lipi.exe <hydrogen bond donor> <hydrogen bond acceptor> <lopP value> <molar refractivity value> <wiener index>\n");
        exit(1);
        }

	float logp=0,re=0, wiener=0;
	int da=0,dd=0;
	char s1[20],s2[20],s3[20],s4[20],s5[20];

	f1=fopen(argv[1],"r");
	f2=fopen(argv[2],"r");
	f3=fopen(argv[3],"r");
	f4=fopen(argv[4],"r");
	f5=fopen(argv[5],"r");

	while(fgets(s1,20,f1)!=NULL)
	dd=atoi(s1);
	while(fgets(s2,20,f2)!=NULL)
	da=atoi(s2);
	while(fgets(s3,20,f3)!=NULL)
	logp=atof(s3);
	while(fgets(s4,20,f4)!=NULL)
	re=atof(s4);
	while(fgets(s5,20,f5)!=NULL)
	wiener=atof(s5);

	printf("%d\t%d\t%.3f\t%.3f\t%.3f\n",dd, da, logp, re, wiener);
	fclose(f1);
	fclose(f2);
	fclose(f3);
	fclose(f4);
	fclose(f5);
	}
