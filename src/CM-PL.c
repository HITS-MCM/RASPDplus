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
 * Purpose: Calculate center of mass of a ligand (drug or drug-like small molecule).
 * Input file: 3D structure of Protein-ligand complex in pdb file format format with identifier ID is always *DRG* without chain informaton
 * How to run:
 * ./CM-PL.exe <PDB-ID> 
 * say,
 * ./pdb2noh.exe ligand.pdb
*/

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include <math.h>
int main(int argc, char *argv[])
{
	FILE *fp;
	if (argc<2) {
  	printf("./CM-PL.exe <PDB-ID> \n");
  	exit(1);
	}
	float sum=0;
	char str[100];

	char atom[6],atmsmb[6],resid[4];
	int atomno,resno;
	float x, y, z, chg, rad, epsi;
	float mx=0,my=0,mz=0;
        float m;
	float com_x=0, com_y=0, com_z=0;

	fp=fopen(argv[1],"r");
	while(fgets(str,100,fp)!=NULL)
	{
		if(strstr(str,"DRG"))
		{
			sscanf(str,"%s %d %s %s %d %f %f %f %f %f %f",atom, &atomno, atmsmb, resid, &resno, &x, &y, &z, &chg, &rad, &epsi);
                        if(strncmp(atmsmb,"C",1)==0 || strncmp(atmsmb,"1C",2)==0 || strncmp(atmsmb,"2C",2)==0 || strncmp(atmsmb,"3C",2)==0)
                                m=12;
                        if(strncmp(atmsmb,"N",1)==0 || strncmp(atmsmb,"1N",2)==0 || strncmp(atmsmb,"2N",2)==0 || strncmp(atmsmb,"3N",2)==0)
                                m=14;
                        if(strncmp(atmsmb,"H",1)==0 || strncmp(atmsmb,"1H",2)==0 || strncmp(atmsmb,"2H",2)==0 || strncmp(atmsmb,"3H",2)==0)
                                m=1;
                        if(strncmp(atmsmb,"O",1)==0 || strncmp(atmsmb,"1O",2)==0 || strncmp(atmsmb,"2O",2)==0 || strncmp(atmsmb,"3O",2)==0)
                                m=16;
                        if(strncmp(atmsmb,"S",1)==0 || strncmp(atmsmb,"1S",2)==0 || strncmp(atmsmb,"2S",2)==0 || strncmp(atmsmb,"3S",2)==0)
                                m=32;
                        if(strncmp(atmsmb,"Se",2)==0 || strncmp(atmsmb,"1Se",3)==0 || strncmp(atmsmb,"2Se",3)==0 || strncmp(atmsmb,"3Se",3)==0)
                                m=79;
                        if(strncmp(atmsmb,"P",1)==0 || strncmp(atmsmb,"1P",2)==0 || strncmp(atmsmb,"2P",2)==0 || strncmp(atmsmb,"3P",2)==0)
                                m=31;
                        if(strncmp(atmsmb,"Cl",2)==0 || strncmp(atmsmb,"CL",2)==0)
                                m=35.5;
                        if(strncmp(atmsmb,"Br",2)==0 || strncmp(atmsmb,"BR",2)==0)
                                m=79;
                        if(strncmp(atmsmb,"F",1)==0)
                                m=19;
                        if(strncmp(atmsmb,"I",1)==0)
                                m=127;


			sum=sum+m;
			mx=mx+m*x;
			my=my+m*y;
			mz=mz+m*z;
		}
	}

com_x=mx/sum;
com_y=my/sum;
com_z=mz/sum;

printf("%8.3f%8.3f%8.3f\n", com_x, com_y, com_z);

fclose(fp);
}
