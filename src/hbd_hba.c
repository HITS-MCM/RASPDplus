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
 * Purpose: Count hydrogen bond donor and acceptor in a molecule. 
 * All hydrogen which are connected with oxygen and nitrogen atoms are donors and all Oxygen and Nitrogen atoms except positively charged and pyrrole-type nitrogen are acceptors.
 * How to run:
 * ./hbd_hba.exe <Output file of lp2.exe> <hydrogen bond donor> <hydrogen bond acceptor> 
*/

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
int main(int argc, char *argv[])
	{
	FILE *f1,*f2,*fr;

	if (argc<2) 
	{
	printf("                   <Input file>          Output file-1           Output file-2\n");	
  	printf("./hbd_hba.exe <Output file of lp2.exe> <hydrogen bond donor> <hydrogen bond acceptor>\n");
  	exit(1);
	}

	char s[80];
	int sum_hbd=0,sum_hba=0;

	fr=fopen(argv[1],"r");
	f1=fopen(argv[2],"w");
	f2=fopen(argv[3],"w");

	if (f1==NULL || f2==NULL)
        {
        printf("hbd_hba.exe: Cannot write output. No output file name is given\n");
        exit(1);
        }

	while(fgets(s,80,fr)!=NULL)
		{
		if(strncmp(s,"HO1",3)==0||strncmp(s,"HN1",3)==0)
		sum_hbd=sum_hbd+1;

		if(strncmp(s,"O",1)==0)
		sum_hba=sum_hba+1;

		if(strncmp(s,"N",1)==0)
			{
			if(strchr(s,'2')==NULL && strchr(s,'4')==NULL)
				{
				if(strlen(s)<=10)
				sum_hba=sum_hba+1;
				}
			}
			if(strncmp(s,"N",1)==0)
			{
			if(strchr(s,'2') || strchr(s,'4'))
				{
				if(strlen(s)<=7)
				sum_hba=sum_hba+1;
				}
			}
		}

fprintf(f1,"%d\n",sum_hbd); /*DONOR*/
fprintf(f2,"%d\n",sum_hba); /*ACCEPTOR*/

fclose(fr);
	}
