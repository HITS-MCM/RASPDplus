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
 * Purpose: This code print predicted binding free energies of small drug-like molecules having physicochemical parameters within the given/default range as supplied by a file name ../data/select_parameter.txt.
 */

#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include<math.h>

int main(int argc, char *argv[])
{
FILE *fp, *fw;

if (argc<15) {
  printf("		 <Inputs>       1         2       3       4        5        6        7         8        9       10       11      12      13    14    <Output>\n");	
  printf("./final_select_all.exe <pred_be_all> <wimin> <wimax> <hbdmin> <hbdmax> <hbamin> <hbamax> <logpmin> <logpmax> <mrmin> <mrmax> <mwmin> <mwmax> <aff> <output file>\n");
  exit(1);
}

float wimin=atof(argv[2]);
float wimax=atof(argv[3]);
int hbdmin=atoi(argv[4]);
int hbdmax=atoi(argv[5]);
int hbamin=atoi(argv[6]);
int hbamax=atoi(argv[7]);
float logpmin=atof(argv[8]);
float logpmax=atof(argv[9]);
float mrmin=atof(argv[10]);
float mrmax=atof(argv[11]);
float mwmin=atoi(argv[12]);
float mwmax=atoi(argv[13]);
float aff=atof(argv[14]);

float W=0,P=0,R=0;
int A=0,D=0;
char s[550],pdbid[20];
char znid[20];
float aff_erf=0, aff_rf=0, aff_dnn=0, aff_knn=0, aff_svr=0, aff_esvr=0, aff_lr=0;
float stdv_erf=0, stdv_rf=0, stdv_dnn=0, stdv_knn=0, stdv_svr=0, stdv_esvr=0, stdv_lr=0;
float mw=0;
fp=fopen(argv[1],"r");
fw=fopen(argv[15],"w");

if (fw==NULL)            
        {
        printf("final_select_all.exe: Cannot write output. No output file name is given\n");
        exit(1);
        }

printf("MoleculeID;PBFE(ERF);PBFE(RF);PBFE(DNN);PBFE(KNN);PBFE(SVR);PBFE(εSVR);PBFE(LR)\n");
	while(fgets(s,500,fp)!=NULL)
	{

	sscanf(s,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %f %f %f %s", &aff_erf, &stdv_erf, &aff_rf, &stdv_rf, &aff_dnn, &stdv_dnn, &aff_knn, &stdv_knn, &aff_svr, &stdv_svr, &aff_esvr, &stdv_esvr, &aff_lr, &stdv_lr, &D, &A, &P, &R, &W, &mw, pdbid);

	   	if(W>wimin && W<wimax    &&
	   	D>=hbdmin && D<=hbdmax   &&
	   	A>=hbamin && A<=hbamax   &&
	   	P>=logpmin && P<=logpmax &&
	   	R>=mrmin && R<=mrmax     &&
	   	mw>=mwmin && mw<=mwmax   &&
	   	aff_erf<=aff && aff_rf<=aff && aff_dnn<=aff && aff_lr<=aff)
		{
//printf("%s;%.2f±%0.2f;%0.2f±%0.2f;%0.2f±%0.2f;%0.2f±%0.2f;%0.2f±%0.2f;%0.2f±%0.2f;%0.2f±%0.2f\n", pdbid, aff_erf, stdv_erf, aff_rf, stdv_rf, aff_dnn, stdv_dnn, aff_knn, stdv_knn, aff_svr, stdv_svr, aff_esvr, stdv_esvr, aff_lr, stdv_lr);
		printf("%s;\t%4.3f; %4.3f; %4.3f; %4.3f; %4.3f; %4.3f; %4.3f\n", pdbid, aff_erf, aff_rf, aff_dnn, aff_knn, aff_svr, aff_esvr, aff_lr);
		fprintf(fw, "%s\n", pdbid);
		}
	}
fclose(fp);
}

