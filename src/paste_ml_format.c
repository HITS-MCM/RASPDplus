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
 * Purpose: This code print the physicochemical parameters seperated by ";" and printing a line with all zero values if of wiener, molar refractivity or, molecular weight column found zero.
 * Input file is the physicochemical parameters of protein and ligands seperated by a space.
 * How to run:
 * ./paste_ml_format.exe Protein-Ligand_Parameter >Protein-Ligand_Parameter.txt
 * Input file (Protein-Ligand_Parameter)
 * 1.538    1.119   0.699   0.839   0.280   0.000   0.000   0.000   1.818   2.098   0.000   0.000   10.468  0.000  1743.19 4       14      -3.519  85.886  1461.920        448.0
 * Output file (Protein-Ligand_Parameter.txt):
 * PDBID;PD(K+R+HIP);PA(D+E);PD(T+S+Y+DH+EH);PA(N+Q+T+S+DH+EH);PD(LYN+N+Q);PA(LYN);PD(W+H);PA(Y+H);PD(Amide-NH);PA(Amide-O);PP(Non-Arom);PP(Arom);PMR(Non-Arom);PMR(Arom);Volume;D;A;P;MR;W;MASS;Expt_BE;Atom_Efficiency;No. of atoms
        Zinc;1.538;1.119;0.699;0.839;0.280;0.000;0.000;0.000;1.818;2.098;0.000;0.000;10.468;0.000;1743.190;4;14;-3.519;85.886;1461.920;448.00;0.000;0.000000;0
 */

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
int main(int argc, char *argv[])
{
FILE *fp;

if (argc<2) {
  printf("./paste_ml_format.exe Protein-Ligand_Parameter\n");
  exit(1);
}

int hbd=0, hba=0, heavy_atom=0;
float logp=0, mr=0, wiener=0, rotbond=0, molW=0, expt_be=0;
char st[500];
char pdbid[5]="Zinc";
float pd_charged=0, pa_charged=0, pd_oh=0, pa_O=0, pd_nh=0, pa_lyn=0, pd_arom=0, pa_arom=0;
float pd_amide=0, pa_amide=0, pp_ali=0, pp_arom=0, mr_ali=0, mr_arom=0, volm=0;
float atom_effi=0;

fp=fopen(argv[1],"r");

printf("PDBID;PD(K+R+HIP);PA(D+E);PD(T+S+Y+DH+EH);PA(N+Q+T+S+DH+EH);PD(LYN+N+Q);PA(LYN);PD(W+H);PA(Y+H);PD(Amide-NH);PA(Amide-O);PP(Non-Arom);PP(Arom);PMR(Non-Arom);PMR(Arom);Volume;D;A;P;MR;W;MASS;Expt_BE;Atom_Efficiency;No. of atoms\n");

while(fgets(st,500,fp)!=NULL)
{
sscanf(st,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %f %f %f", &pd_charged, &pa_charged, &pd_oh, &pa_O, &pd_nh, &pa_lyn, &pd_arom, &pa_arom, &pd_amide, &pa_amide, &pp_ali, &pp_arom, &mr_ali, &mr_arom, &volm, &hbd, &hba, &logp, &mr, &wiener, &molW);

//if(hbd==0 && hba==0 && logp==0 && mr==0 && wiener==0 && molW==0)
if((logp==0 && mr==0) || wiener==0 || molW==0)
printf("0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0\n");
else
printf("%s;%0.3f;%0.3f;%0.3f;%0.3f;%0.3f;%0.3f;%0.3f;%0.3f;%0.3f;%0.3f;%0.3f;%0.3f;%0.3f;%0.3f;%0.3f;%d;%d;%0.3f;%0.3f;%0.3f;%0.2f;%0.3f;%0.6f;%d\n", pdbid, pd_charged, pa_charged, pd_oh, pa_O, pd_nh, pa_lyn, pd_arom, pa_arom, pd_amide, pa_amide, pp_ali, pp_arom, mr_ali, mr_arom, volm, hbd, hba, logp, mr, wiener, molW, expt_be, atom_effi, heavy_atom);
}
fclose(fp);
}
