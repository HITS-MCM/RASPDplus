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
 * Purpose: This code print the distance of an atom in a molecule from its center of mass.
 * How to run:
 * ./maxD_calculator.exe <ligand id> <com_ligand> | sort -n | tail -n 1 
 * ./maxD_calculator.exe ligand.pdb com_ligand | | sort -n | tail -n 1 [say center of mass of the ligand is 0,0,0 here]
 */

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include <math.h>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
using namespace std;
typedef struct protein {    
        float x,y,z;
        int resno, atomno;
        char  atomname[5];
        char resname[4];
	char atom[6];
}P;
typedef struct COMass {  
        float x,y,z;
}COM;


int main(int argc, char *argv[])
{
        FILE *fp,*fd;
	if (argc<3) {
	printf("                   <Input file-1> <Input file-2>\n");
  	printf("./maxD_calculator.exe ligand.pdb com_ligand \n");
  	exit(1);
	}

        char s[100];
	double distance=0;
        fp=fopen(argv[1],"r");
        fd=fopen(argv[2],"r");
        int counter_protein=0, counter_com=0;
        P   P1 [70000];
        COM P2 [70000];
        vector <string> PDB;

        while (fgets (s, 100, fp)!=NULL){
                if (8==sscanf(s,"%s%d%s%s%d%f%f%f",P1[counter_protein].atom, &P1[counter_protein].atomno, P1[counter_protein].atomname\
                                            ,P1[counter_protein].resname,&P1[counter_protein].resno,&P1[counter_protein].x\
                                            ,&P1[counter_protein].y,&P1[counter_protein].z)){
                counter_protein++;
                }
        }
        while (fgets (s, 100, fd)!=NULL){
                if(3==sscanf(s,"%f%f%f",&P2[counter_com].x, &P2[counter_com].y, &P2[counter_com].z)){
		counter_com++;
                }
        }
        int i=0, j=0;
        for (i=0;i<counter_protein;i++)
		{
                		for ( j=0;j<counter_com;j++)
				{
                        		distance=sqrt((P1[i].x - P2[j].x) * (P1[i].x - P2[j].x) + (P1[i].y - P2[j].y) * (P1[i].y - P2[j].y) + (P1[i].z - P2[j].z) * (P1[i].z - P2[j].z) );
					{
					printf("%0.3lf\n", distance);
					}
        		        }
        
		}

}
