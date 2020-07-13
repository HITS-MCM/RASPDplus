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
 * Purpose: This program differentiate atom types of different types.
 * Say, in Toluene (a benzene ring with one methyl group attach, Ph-CH3), all the six carbons are not similar. 5 Carbons of benzne ring are same and the carbon of benzene where, methyl group is attach is different than the other 5 carbons. Hence, to differentiate this carbon atom where methyl group is attached, this program add a new three letter nomenclature, "CMe", with the benzene carbon.
 In order to differentiate the methyl group which is attach with non-aromatic and aromatic group, this program will add a new three letter nomenclature, "CC4", with the benzylic carbon.
 Thus, this program will differentiate atom types of different types.
 * How to run:
 * ./lp5.exe <bond order of each of the atoms in a molecule, outout of lp4.pl program> <Output file of lp2.exe> <Output file of lp3.exe> <Output file name> 
*/

#include<stdio.h>
#include<string.h>
#include<stdlib.h>

int main(int argc, char *argv[])
{
	FILE *fp,*fr,*fw,*f1;
	  if (argc<4)
          {
          printf("                                <Input file-1>                                          <Input file-2>           <Input file-3>           <Output file>\n");
          printf("./lp5.exe <bond order of each of the atoms in a molecule, outout of lp4.pl program> <Output file of lp2.exe> <Output file of lp3.exe> <Output file name>\n");
          exit(1);
          }
       	struct format
       {
        int a;
        int b;
        int c;
	int d;
	int e;
	int f;
	int g;	
       };
        struct format f;

        char strings[500],str[500],s[500],st[50];
	char strings_array[700][50];
	char str_array[500][50],bnd[700][50];
	char final_string[500][50];
 	char *result;	
	int bo=0,ln=0,ln1,ln2,ln3,l;
	f1=fopen(argv[1],"r");
	fp=fopen(argv[2],"r");
	fr=fopen(argv[3],"r");
	fw=fopen(argv[4],"w");
	if (fw==NULL)            
        {
        printf("lp5.exe: Cannot write output. No output file name is given\n");
        exit(1);
        }	
	//int ln=0;
	while(fgets(st,50,f1)!=NULL)
	{
	strncpy(bnd[ln],st,strlen(st)-1);
	ln++;
	}
	
	int line1=0;
	while(fgets(strings,500,fp)!=NULL)
	{
		int len = strlen(strings);
		strncpy(strings_array[line1],strings,len-1);
		//printf("%s\n",strings_array[line1]);
		line1++;
	}

	int line=0;
	while(fgets(str,500,fr)!=NULL)
	{
                f.a=0;
                f.b=0;
                f.c=0;
		f.d=0;
		f.e=0;
		f.f=0;
		f.g=0;
                sscanf(str,"%d%d%d%d%d%d%d",&f.a,&f.b,&f.c,&f.d,&f.e,&f.f,&f.g);
                l=strlen(strings_array[line]);
		bo=atoi(bnd[line]);
		if(strncmp("C",strings_array[f.a-1],1)==0)
		{
		if(l==12){
		if((strchr(strings_array[f.b-1],'4') && strlen(strings_array[f.b-1])==9)||(strchr(strings_array[f.c-1],'4') && strlen(strings_array[f.c-1])==9)||(strchr(strings_array[f.d-1],'4') && strlen(strings_array[f.d-1])==9)||(strchr(strings_array[f.e-1],'4') && strlen(strings_array[f.e-1])==9))
		strcat(strings_array[f.a-1],"CC4");
		if((strchr(strings_array[f.b-1],'4') && strlen(strings_array[f.b-1])==12 && (strstr(strings_array[f.b-1],"CX4")||strstr(strings_array[f.b-1],"CMe")))||(strchr(strings_array[f.c-1],'4') && strlen(strings_array[f.c-1])==12 && (strstr(strings_array[f.c-1],"CX4")||strstr(strings_array[f.c-1],"CMe")))||(strchr(strings_array[f.d-1],'4') && strlen(strings_array[f.d-1])==12 && (strstr(strings_array[f.d-1],"CX4")||strstr(strings_array[f.d-1],"CMe")))||(strchr(strings_array[f.e-1],'4') && strlen(strings_array[f.e-1])==12 && (strstr(strings_array[f.e-1],"CX4")||strstr(strings_array[f.e-1],"CMe")))){
		if(strstr(strings_array[f.a-1],"CC4")==NULL)
		strcat(strings_array[f.a-1],"CC4");}}


		if(l==9 && strchr(strings_array[f.a-1],'4')==NULL){
		//if(strchr(strings_array[f.b-1],'4')||strchr(strings_array[f.c-1],'4')||strchr(strings_array[f.d-1],'4')||strchr(strings_array[f.e-1],'4'))
		if((strchr(strings_array[f.b-1],'4') && strlen(strings_array[f.b-1])==9)||(strchr(strings_array[f.c-1],'4') && strlen(strings_array[f.c-1])==9)||(strchr(strings_array[f.d-1],'4') && strlen(strings_array[f.d-1])==9)||(strchr(strings_array[f.e-1],'4') && strlen(strings_array[f.e-1])==9))
		strcat(strings_array[f.a-1],"CC4");
		if((strchr(strings_array[f.b-1],'4') && strlen(strings_array[f.b-1])==12 && (strstr(strings_array[f.b-1],"CMe")||strstr(strings_array[f.b-1],"CX4")))||(strchr(strings_array[f.c-1],'4') && strlen(strings_array[f.c-1])==12 &&(strstr(strings_array[f.c-1],"CMe")||strstr(strings_array[f.c-1],"CX4")))||(strchr(strings_array[f.d-1],'4') && strlen(strings_array[f.d-1])==12 && (strstr(strings_array[f.d-1],"CMe")||strstr(strings_array[f.d-1],"CX4")))||(strchr(strings_array[f.e-1],'4') && strlen(strings_array[f.e-1])==12 && (strstr(strings_array[f.e-1],"CMe")||strstr(strings_array[f.e-1],"CX4")))){
		if(strstr(strings_array[f.a-1],"CC4")==NULL)
		strcat(strings_array[f.a-1],"CC4");}}

		if(l==9 && strchr(strings_array[f.a-1],'4') && strstr(strings_array[f.a-1],"CC1") && strchr(strings_array[f.a-1],'2')==NULL){
		ln1=strlen(strings_array[f.b-1]);
		ln2=strlen(strings_array[f.c-1]);
		ln3=strlen(strings_array[f.d-1]);
		if(ln1==9 && ln2==9 && ln3==9 && strchr(strings_array[f.b-1],'2')==NULL && strchr(strings_array[f.c-1],'2')==NULL && strchr(strings_array[f.d-1],'2')==NULL && strchr(strings_array[f.b-1],'3')==NULL && strchr(strings_array[f.c-1],'3')==NULL && strchr(strings_array[f.d-1],'3')==NULL)
		strcat(strings_array[f.a-1],"CX4");

		if(((ln1==12&&strchr(strings_array[f.b-1],'2')==NULL&&strchr(strings_array[f.b-1],'4')&&(strstr(strings_array[f.b-1],"CMe")||strstr(strings_array[f.b-1],"CX4"))) || (ln2==12&&strchr(strings_array[f.c-1],'2')==NULL&&strchr(strings_array[f.c-1],'4')&&(strstr(strings_array[f.c-1],"CMe")||strstr(strings_array[f.c-1],"CX4"))) || (ln3==12&&strchr(strings_array[f.d-1],'2')==NULL&&strchr(strings_array[f.d-1],'4')&&(strstr(strings_array[f.d-1],"CMe")||strstr(strings_array[f.d-1],"CX4"))))&&strchr(strings_array[f.b-1],'4')&&strchr(strings_array[f.c-1],'4')&&strchr(strings_array[f.d-1],'4')){
		if(strstr(strings_array[f.a-1],"CX4")==NULL)
		strcat(strings_array[f.a-1],"CX4");}

		if((ln1>=12&&strchr(strings_array[f.b-1],'2')==NULL&&strstr(strings_array[f.b-1],"CMe")==NULL&&strstr(strings_array[f.b-1],"CX4")==NULL)||(ln2>=12&&strchr(strings_array[f.c-1],'2')==NULL&&strstr(strings_array[f.c-1],"CMe")==NULL&&strstr(strings_array[f.c-1],"CX4")==NULL)||(ln3>=12&&strchr(strings_array[f.d-1],'2')==NULL&&strstr(strings_array[f.d-1],"CMe")==NULL&&strstr(strings_array[f.d-1],"CX4")==NULL))
		strcat(strings_array[f.a-1],"CMe");
		if((ln1==12&&strchr(strings_array[f.b-1],'2')&&strstr(strings_array[f.b-1],"CC4")&&strstr(strings_array[f.b-1],"CMe")==NULL&&strstr(strings_array[f.b-1],"CX4")==NULL)||(ln2==12&&strstr(strings_array[f.c-1],"CC4")&&strchr(strings_array[f.c-1],'2')&&strstr(strings_array[f.c-1],"CMe")==NULL&&strstr(strings_array[f.c-1],"CX4")==NULL)||(ln3==12&&strstr(strings_array[f.d-1],"CC4")&&strchr(strings_array[f.d-1],'2')&&strstr(strings_array[f.d-1],"CMe")==NULL&&strstr(strings_array[f.d-1],"CX4")==NULL)){
		if(strstr(strings_array[f.a-1],"CMe")==NULL)
		strcat(strings_array[f.a-1],"CMe");}
		
		if(((ln1==9&&strchr(strings_array[f.b-1],'2')) ||( ln2==9&&strchr(strings_array[f.c-1],'2')) || (ln3==9&&strchr(strings_array[f.d-1],'2'))) && strstr(strings_array[f.a-1],"CMe")==NULL)
		strcat(strings_array[f.a-1],"CMe");

		if((ln1==9&&strchr(strings_array[f.b-1],'3')&&strstr(strings_array[f.b-1],"CC4")&&strstr(strings_array[f.b-1],"CMe")==NULL&&strstr(strings_array[f.b-1],"CX4")==NULL)||(ln2==9&&strstr(strings_array[f.c-1],"CC4")&&strchr(strings_array[f.c-1],'3')&&strstr(strings_array[f.c-1],"CMe")==NULL&&strstr(strings_array[f.c-1],"CX4")==NULL)||(ln3==9&&strstr(strings_array[f.d-1],"CC4")&&strchr(strings_array[f.d-1],'3')&&strstr(strings_array[f.d-1],"CMe")==NULL&&strstr(strings_array[f.d-1],"CX4")==NULL)){
                if(strstr(strings_array[f.a-1],"CMe")==NULL)
                strcat(strings_array[f.a-1],"CMe");}

		if(((ln1==6 &&strchr(strings_array[f.b-1],'3'))|| (ln2==6&&strchr(strings_array[f.c-1],'3')) || (ln3==6&&strchr(strings_array[f.d-1],'3'))) && strstr(strings_array[f.a-1],"CMe")==NULL)
		strcat(strings_array[f.a-1],"CMe");}
		//if(strstr(strings_array[f.a-1],"CO2") || strstr(strings_array[f.a-1],"CN2") || strstr(strings_array[f.a-1],"CS2") || strstr(strings_array[f.a-1],"CP2")){
		//if(strchr(strings_array[f.b-1],'4')||strchr(strings_array[f.c-1],'4')||strchr(strings_array[f.d-1],'4')||strchr(strings_array[f.e-1],'4'))
		//strcat(strings_array[f.a-1],"CC4");}
		}

		if(strncmp("N",strings_array[f.a-1],1)==0)
                {
                if(strchr(strings_array[f.a-1],'4')==NULL && strchr(strings_array[f.a-1],'2')==NULL && strchr(strings_array[f.a-1],'3')==NULL)
		{
		if(l==12)
		strcat(strings_array[f.a-1],"PVE");
                }
		if(strchr(strings_array[f.a-1],'2') && strchr(strings_array[f.a-1],'4')==NULL && strchr(strings_array[f.a-1],'3')==NULL)
		{
		if(l==9)
		strcat(strings_array[f.a-1],"PVE");
                }
		if(strchr(strings_array[f.a-1],'4') && l==9 && bo==4)
		strcat(strings_array[f.a-1],"PVE");
		
		if(l==6 &&(strchr(strings_array[f.a-1],'4')==NULL && strchr(strings_array[f.a-1],'2')==NULL && strchr(strings_array[f.a-1],'3')==NULL))
		strcat(strings_array[f.a-1],"NVE");
		if(l==3 && strchr(strings_array[f.a-1],'2'))
		strcat(strings_array[f.a-1],"NVE");
		if(l==6 && strchr(strings_array[f.a-1],'4') && bo==2)
		strcat(strings_array[f.a-1],"NVE");
		
		if(strchr(strings_array[f.a-1],'4')==NULL){
                if((strstr(strings_array[f.b-1],"CC4")&&strncmp("C",strings_array[f.b-1],1)==0&&strlen(strings_array[f.b-1])==9)||(strstr(strings_array[f.c-1],"CC4")&&strncmp("C",strings_array[f.c-1],1)==0&&strlen(strings_array[f.c-1])==9)||(strstr(strings_array[f.d-1],"CC4")&&strncmp("C",strings_array[f.d-1],1)==0&&strlen(strings_array[f.d-1])==9)||strstr(strings_array[f.b-1],"CN4")||strstr(strings_array[f.c-1],"CN4")||strstr(strings_array[f.d-1],"CN4")||strstr(strings_array[f.b-1],"CO4")||strstr(strings_array[f.c-1],"CO4")||strstr(strings_array[f.d-1],"CO4")||strstr(strings_array[f.b-1],"CS4")||strstr(strings_array[f.c-1],"CS4")||strstr(strings_array[f.d-1],"CS4")||strstr(strings_array[f.b-1],"NN4")||strstr(strings_array[f.c-1],"NN4")||strstr(strings_array[f.d-1],"NN4")||strstr(strings_array[f.b-1],"NS4")||strstr(strings_array[f.c-1],"NS4")||strstr(strings_array[f.d-1],"NS4")||strstr(strings_array[f.b-1],"NC4")||strstr(strings_array[f.c-1],"NC4")||strstr(strings_array[f.d-1],"NC4")||strstr(strings_array[f.b-1],"NO4")||strstr(strings_array[f.c-1],"NO4")||strstr(strings_array[f.d-1],"NO4"))
		strcat(strings_array[f.a-1],"CC4");}
		
		}
		//////////////////////////////////////////
		if(strncmp(strings_array[f.a-1],"O",1)==0)
		{
		if(strstr(strings_array[f.a-1],"OH1") || strstr(strings_array[f.a-1],"OC1") || strstr(strings_array[f.a-1],"ON1") || strstr(strings_array[f.a-1],"OO1") || strstr(strings_array[f.a-1],"OF1") || strstr(strings_array[f.a-1],"OP1") || strstr(strings_array[f.a-1],"OS1") || strstr(strings_array[f.a-1],"Ol1") || strstr(strings_array[f.a-1],"OB1") || strstr(strings_array[f.a-1],"OI1")){  //for carbonyl group attach neutral acidic oxygen.
                if(l==6){
                if(strstr(strings_array[f.b-1],"CO2") || strstr(strings_array[f.c-1],"CO2") || strstr(strings_array[f.b-1],"CN2") || strstr(strings_array[f.c-1],"CN2") || strstr(strings_array[f.b-1],"CS2") || strstr(strings_array[f.c-1],"CS2") || strstr(strings_array[f.b-1],"CC2") || strstr(strings_array[f.c-1],"CC2"))
                strcat(strings_array[f.a-1],"CO2");}}

		if(strstr(strings_array[f.a-1],"OH1")==NULL && strstr(strings_array[f.a-1],"OC4")==NULL && strstr(strings_array[f.a-1],"ON4")==NULL && strstr(strings_array[f.a-1],"OS4")==NULL && strstr(strings_array[f.a-1],"OO4")==NULL && strstr(strings_array[f.a-1],"OP4")==NULL){  //for aromatic etheric oxygen.
		//if(strstr(strings_array[f.b-1],"O2")==NULL && strstr(strings_array[f.c-1],"O2")==NULL && strstr(strings_array[f.b-1],"N2")==NULL && strstr(strings_array[f.c-1],"N2")==NULL && strstr(strings_array[f.b-1],"S2")==NULL && strstr(strings_array[f.c-1],"S2")==NULL){ //this is nullified because esteric group is consider as etheric type.
		//if(strchr(strings_array[f.b-1],'4') || strchr(strings_array[f.c-1],'4'))
		if((strstr(strings_array[f.b-1],"CC4")&&strncmp("C",strings_array[f.b-1],1)==0&&strlen(strings_array[f.b-1])==9)||(strstr(strings_array[f.c-1],"CC4")&&strncmp("C",strings_array[f.c-1],1)==0&&strlen(strings_array[f.c-1])==9)||(strstr(strings_array[f.d-1],"CC4")&&strncmp("C",strings_array[f.d-1],1)==0&&strlen(strings_array[f.d-1])==9)||strstr(strings_array[f.b-1],"CN4")||strstr(strings_array[f.c-1],"CN4")||strstr(strings_array[f.d-1],"CN4")||strstr(strings_array[f.b-1],"CO4")||strstr(strings_array[f.c-1],"CO4")||strstr(strings_array[f.d-1],"CO4")||strstr(strings_array[f.b-1],"CS4")||strstr(strings_array[f.c-1],"CS4")||strstr(strings_array[f.d-1],"CS4")||strstr(strings_array[f.b-1],"NN4")||strstr(strings_array[f.c-1],"NN4")||strstr(strings_array[f.d-1],"NN4")||strstr(strings_array[f.b-1],"NS4")||strstr(strings_array[f.c-1],"NS4")||strstr(strings_array[f.d-1],"NS4")||strstr(strings_array[f.b-1],"NC4")||strstr(strings_array[f.c-1],"NC4")||strstr(strings_array[f.d-1],"NC4")||strstr(strings_array[f.b-1],"NO4")||strstr(strings_array[f.c-1],"NO4")||strstr(strings_array[f.d-1],"NO4"))
		strcat(strings_array[f.a-1],"CC4");}

		if(strstr(strings_array[f.a-1],"ON1")){  //for nitrate oxygen.
		if(strstr(strings_array[f.b-1],"NO2")||strstr(strings_array[f.c-1],"NO2"))
		strcat(strings_array[f.a-1],"NO2");}
		
		if(strstr(strings_array[f.a-1],"OC2") && strchr(strings_array[f.b-1],'4') && strlen(strings_array[f.b-1])==9)  //for aromatic carbonyl oxygen.
		strcat(strings_array[f.a-1],"Oc4");
		
		if(strstr(strings_array[f.a-1],"OC2")){  //for aliphatic/aromatic carboxylic oxygen.
		if(strstr(strings_array[f.b-1],"CC1CC1") || strstr(strings_array[f.b-1],"CC1CH1")||strstr(strings_array[f.b-1],"CC1CO2CH1")||strstr(strings_array[f.b-1],"CH1CC1") || strstr(strings_array[f.b-1],"CC1CN1")||strstr(strings_array[f.b-1],"CN1CC1") || strstr(strings_array[f.b-1],"CC1CO1") ||strstr(strings_array[f.b-1],"CC1CO2CO1") ||strstr(strings_array[f.b-1],"CO1CC1") || strstr(strings_array[f.b-1],"CC1CS1")||strstr(strings_array[f.b-1],"CC1CO2CS1")||strstr(strings_array[f.b-1],"CS1CC1") || strstr(strings_array[f.b-1],"CC1CF1")||strstr(strings_array[f.b-1],"CC1CO2CF1")||strstr(strings_array[f.b-1],"CF1CC1") || strstr(strings_array[f.b-1],"CC1Cl1")||strstr(strings_array[f.b-1],"CC1CO2Cl1")||strstr(strings_array[f.b-1],"Cl1CC1") || strstr(strings_array[f.b-1],"CC1CB1")||strstr(strings_array[f.b-1],"CC1CO2CB1")||strstr(strings_array[f.b-1],"CB1CC1") || strstr(strings_array[f.b-1],"CC1CI1")||strstr(strings_array[f.b-1],"CC1CO2CI1")||strstr(strings_array[f.b-1],"CI1CC1") || strstr(strings_array[f.b-1],"CH1CN1")||strstr(strings_array[f.b-1],"CN1CH1") ||strstr(strings_array[f.b-1],"CN1CO2CH1")|| strstr(strings_array[f.b-1],"CH1CO1")||strstr(strings_array[f.b-1],"CO1CH1") || strstr(strings_array[f.b-1],"CO1CO2CH1")|| strstr(strings_array[f.b-1],"CH1CH1") || strstr(strings_array[f.b-1],"CO2CO2")){
		if(strchr(strings_array[f.b-1],'4')==NULL) //for atleast one carbon and other heteroatom attach aliphatic carboxylic oxygen.
		strcat(strings_array[f.a-1],"CAL");}
		
		if(strstr(strings_array[f.b-1],"CC1CC1") || strstr(strings_array[f.b-1],"CC1CH1")||strstr(strings_array[f.b-1],"CC1CO2CH1")||strstr(strings_array[f.b-1],"CH1CC1") || strstr(strings_array[f.b-1],"CC1CN1")||strstr(strings_array[f.b-1],"CC1CO2CN1")||strstr(strings_array[f.b-1],"CN1CC1") || strstr(strings_array[f.b-1],"CC1CO1")||strstr(strings_array[f.b-1],"CC1CO2CO1")||strstr(strings_array[f.b-1],"CO1CC1") || strstr(strings_array[f.b-1],"CC1CS1")||strstr(strings_array[f.b-1],"CC1CO2CS1")||strstr(strings_array[f.b-1],"CS1CC1") || strstr(strings_array[f.b-1],"CC1CF1")||strstr(strings_array[f.b-1],"CC1CO2CF1")||strstr(strings_array[f.b-1],"CF1CC1") || strstr(strings_array[f.b-1],"CC1Cl1")||strstr(strings_array[f.b-1],"CC1CO2Cl1")||strstr(strings_array[f.b-1],"Cl1CC1") || strstr(strings_array[f.b-1],"CC1CB1")||strstr(strings_array[f.b-1],"CC1CO2CB1")||strstr(strings_array[f.b-1],"CB1CC1") || strstr(strings_array[f.b-1],"CC1CI1")||strstr(strings_array[f.b-1],"CC1CO2CI1")||strstr(strings_array[f.b-1],"CI1CC1")){  ////for atleast one aromatic carbon and other heteroatom attach carboxylic oxygen.
		if(strchr(strings_array[f.b-1],'4'))  
                strcat(strings_array[f.a-1],"CAR");}  
		if(strstr(strings_array[f.b-1],"CC1")==NULL && strstr(strings_array[f.b-1],"CH1")==NULL && (strstr(strings_array[f.b-1],"CN1") || strstr(strings_array[f.b-1],"CO1") || strstr(strings_array[f.b-1],"CP1") || strstr(strings_array[f.b-1],"CS1") || strstr(strings_array[f.b-1],"Cl1") || strstr(strings_array[f.b-1],"CF1") || strstr(strings_array[f.b-1],"CB1") || strstr(strings_array[f.b-1],"CI1"))) //for both heteroatom attach aliphatic/aromatic carboxylic oxygen.
		strcat(strings_array[f.a-1],"CHE");}
		
		if(l==3 && strncmp(strings_array[f.a-1],"OC1",3)==0){  //For carboxylate oxygen only.
		if(strstr(strings_array[f.b-1],"CO2"))
		strcat(strings_array[f.a-1],"CO2");}

		//if(strstr(strings_array[f.a-1],"OP1") || strstr(strings_array[f.a-1],"OS1") || strstr(strings_array[f.a-1],"Ol1") || strstr(strings_array[f.a-1],"OB1") || strstr(strings_array[f.a-1],"OI1")){ //for P=O/S=O attach charged acidic oxygen.
                //if(l==3){
		//if(strstr(strings_array[f.b-1],"PO2") || strings_array[f.c-1],"PO2" || strstr(strings_array[f.b-1],"SO2") || strings_array[f.c-1],"SO2")
		//strcat(strings_array[f.a-1],"CO2");}}
		}
		///////////////////////////////////////////
		if(strncmp(strings_array[f.a-1],"H",1)==0)
		{
		if(strstr(strings_array[f.a-1],"HO1")){
		if(strstr(strings_array[f.b-1],"CO2") || strstr(strings_array[f.b-1],"OO1")||strstr(strings_array[f.b-1],"OS1"))
		strcat(strings_array[f.a-1],"ACD");}
		if((strstr(strings_array[f.b-1],"OC1") || strstr(strings_array[f.b-1],"OP1")) && strstr(strings_array[f.b-1],"ON1")==NULL && strstr(strings_array[f.b-1],"OO1")==NULL && strstr(strings_array[f.b-1],"OS1")==NULL && strstr(strings_array[f.b-1],"CO2")==NULL)
                strcat(strings_array[f.a-1],"OOL");
		if(strstr(strings_array[f.a-1],"HN1")){
		if(strstr(strings_array[f.b-1],"ON1"))	
		strcat(strings_array[f.a-1],"AMN");}
		}
		
		strcpy(final_string[line],strings_array[line]);
		
		fprintf(fw,"%s\n",final_string[line]);
                line++;
                }

		     
		fclose(fp);
		fclose(fr);
		fclose(fw);

}


