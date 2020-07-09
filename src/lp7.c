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
 * Purpose: This program will print the atom types of atoms in a molecule as per nomenclature given by Scott A. Wildman and Gordon M. Crippen in their work (J. Chem. Inf. Comput. Sci. 1999, 39, 868).
 * How to run:
 * ./lp7.exe <Output file of lp6.exe> <Output file name> 
*/

#include<stdio.h>
#include<string.h>
#include<stdlib.h>

int main(int argc, char *argv[])
{
	FILE *fr,*fw;
	if (argc<2) {
	printf("               <Input file-1>         <Output file>\n");
  	printf("./lp7.exe <Output file of lp6.exe> <Output file name>\n");
  	exit(1);
	}

        char strings[500],str[500],s[500];
	char strings_array[500][500];
	char str_array[500][500];
	char final_string[500][500];
 	char *result;
	int k,sum=0,sum2=0;
		
	fr=fopen(argv[1],"r");
	fw=fopen(argv[2],"w");
	if (fw==NULL)            
        {
        printf("lp7.exe: Cannot write output. No output file name is given\n");
        exit(1);
        }
	int line=0;
	while(fgets(strings,500,fr)!=NULL)
	{
	sscanf(strings,"%s",s);
        k=strlen(s)+1;
	
	if(strncmp(s,"N",1)==0)
	{
	if(strstr(s,"PVE")==NULL && strchr(s,'2')==NULL){
	if(k==10 && strstr(s,"NH1NH1") && strchr(s,'4')==NULL){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tN1  ");
	printf("N1\n");}
	if(k==13 && strstr(s,"NH1NH1") && strstr(s,"CC4")){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tN3  ");
	printf("N3\n");}
	
        if(k==10 && strstr(s,"NH1") && strstr(s,"NH1NH1")==NULL && strchr(s,'4')==NULL){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tN2  ");
        printf("N2\n");}
        if(k==13 && strstr(s,"NH1") && strstr(s,"NH1NH1")==NULL && strstr(s,"CC4")){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tN4  ");
        printf("N4\n");}

	if(k==10 && strstr(s,"NH1")==NULL && strchr(s,'4')==NULL){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tN7  ");
	printf("N7\n");}
	if(k==13 && strstr(s,"NH1")==NULL && strstr(s,"CC4")){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tN8  ");
	printf("N8\n");}}
	
	if((k==7 || k==10) && (strstr(s,"NH1") && strstr(s,"PVE")==NULL && strchr(s,'2') && strstr(s,"N2  ")==NULL)){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tN5  ");
        printf("N5\n");}
        if((k==7 || k==10) && (strstr(s,"NH1")==NULL && strstr(s,"PVE")==NULL && strchr(s,'2') && strstr(s,"N2  ")==NULL)){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tN6  ");
        printf("N6\n");}

	if(k==4 && strchr(s,'3')){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tN9  ");
	printf("N9\n");}

	if(k>=16 && strstr(s,"PVE") && strstr(s,"NH1") && strstr(s,"NC4")==NULL && strstr(s,"NN4")==NULL && strstr(s,"NO4")==NULL && strstr(s,"NS4")==NULL){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\tN10  ");
	printf("N10\n");}
	if(k==13 && strchr(s,'2') && strstr(s,"NH1") && strstr(s,"PVE")){      //Newly added, not in at_des.c programme.
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tN10  ");
        printf("N10\n");}
	
	if(strstr(s,"NC4")||strstr(s,"NN4")||strstr(s,"NO4")||strstr(s,"NS4")||strstr(s,"NP4")){
	if(strstr(s,"PVE")==NULL){
	if(k==7 || k==10){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tN11  ");
	printf("N11\n");}}}
	
	if(strstr(s,"NC4")||strstr(s,"NN4")||strstr(s,"NO4")||strstr(s,"NS4")||strstr(s,"NP4")){
	if(strstr(s,"PVE")){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tN12  ");	
	printf("N12\n");}}
	
	if(k==16 && strstr(s,"PVE") && strstr(s,"NH1")==NULL && strchr(s,'4')==NULL && strchr(s,'2')==NULL) {
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\tN13  ");
	printf("N13\n");}
	if(k==16 && strstr(s,"PVE") && strstr(s,"NH1")==NULL && strstr(s,"CC4") && strchr(s,'2')) {   //Newly added, not in at_des.c programme.
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\tN13  ");
        printf("N13\n");}
        if(k==13 && strstr(s,"PVE") && strstr(s,"NH1")==NULL && strchr(s,'2') && strchr(s,'4')==NULL){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tN13  ");
        printf("N13\n");}
	if(strstr(s,"NC2NN2") || strstr(s,"NN2NC2") || strstr(s,"NN3NC2")){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tN13  ");
        printf("N13\n");}
	
	if(k==7 && strchr(s,'3') && strchr(s,'2')==NULL && strchr(s,'4')==NULL){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tN14  ");
	printf("N14\n");}
	if(strstr(s,"NVE")){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tN14  ");
        printf("N14\n");}
	if(strstr(s,"NN2NN2") || strstr(s,"NN3NN2") || strstr(s,"NN3NN1")){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tN14  ");
        printf("N14\n");}
	}

	/***********************************************************************/		
	if(strncmp(s,"O",1)==0)
	{
	if(strstr(s,"OC4")||strstr(s,"ON4")||strstr(s,"OS4")||strstr(s,"OP4")){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tO1  ");	
	printf("O1\n");}
	
	if(strstr(s,"OH1") && k>=7){  //for alcoholic oxygen including carbonyl group attach neutral oxygen.
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tO2  ");
	printf("O2\n");}
	if(strstr(s,"OH1OH1")){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tO2  ");
        printf("O2\n");}
	
	if(k==7 && strstr(s,"OH1")==NULL && strchr(s,'4')==NULL && strchr(s,'2')==NULL){   //for ether type including carbonyl group attached oxygen.
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tO3  ");
        printf("O3\n");}
	if(k==10 && strstr(s,"OH1")==NULL && strchr(s,'4')==NULL && strstr(s,"CO2")){   //for ether type including carbonyl group attached oxygen.
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tO3  ");
        printf("O3\n");}
	if(k==10 && strstr(s,"OH1")==NULL  && strchr(s,'4') && strchr(s,'2')==NULL){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tO4  ");
	printf("O4\n");}
	if(k==13 && strstr(s,"OH1")==NULL  && strchr(s,'4') && strstr(s,"CO2")){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tO4  ");
        printf("O4\n");}
	
	if(strstr(s,"OO2")||strstr(s,"ON2")||strstr(s,"ON1NO2")){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tO5  ");
	printf("O5\n");}
	
	if(k==4 && strstr(s,"OS1") ){         //change on 28/08/2012
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tO6  ");
	printf("O6\n");}

	if(k==4 && strstr(s,"OP1") ){        //change on 28/08/2012
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tO7  ");
	printf("O7\n");}
	if(k==4 && strstr(s,"OC1")){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tO7  ");
	printf("O7\n");}
	if(k==4 && (strstr(s,"OF1")||strstr(s,"Ol1")||strstr(s,"OB1")||strstr(s,"OI1")||strstr(s,"OH1")||strstr(s,"OO1"))){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tO7  ");
	printf("O7\n");}
	
	if(strstr(s,"Oc4")){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tO8  ");
	printf("O8\n");}
	
	if(strstr(s,"OC2") && strstr(s,"CAL")){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tO9  ");
	printf("O9\n");}
	
	if(strstr(s,"OC2") && strstr(s,"CAR") && k==7){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tO10  ");
	printf("O10\n");}
        
	if(strstr(s,"OC2") && strstr(s,"CHE")&& k==7){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tO11  ");
	printf("O11\n");}
	
	if(k==7 && strstr(s,"OC1CO2")){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tO12  ");
	printf("O12\n");}
	}
	/***************************************/
	if(strncmp(s,"C",1)==0)
	{
	if(strstr(s,"CH1CH1CH1CH1") && k==13){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC1  ");
	printf("C1\n");}
	if(strstr(s,"CH1CH1CH1") && strstr(s,"CC1") && k==13){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC1  ");
	printf("C1\n");}
	if(strstr(s,"CH1CH1") && strstr(s,"CC1CC1") && k==13){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC1  ");
	printf("C1\n");}

	if(strstr(s,"CH1") && strstr(s,"CC1CC1CC1") && k==13){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC2  ");
	printf("C2\n");}
	if(strstr(s,"CC1CC1CC1CC1") && k==13){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC2  ");
	printf("C2\n");}
	
	if(strstr(s,"CH1CH1CH1") && strstr(s,"CH1CH1CH1CH1")==NULL && k==13){
	if(strstr(s,"CN1")||strstr(s,"CO1")||strstr(s,"CP1")||strstr(s,"CS1")||strstr(s,"Cl1")||strstr(s,"CF1")||strstr(s,"CB1")||strstr(s,"CI1")){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC3  ");
	printf("C3\n");}}
	if(strstr(s,"CH1CH1") && strstr(s,"CH1CH1CH1CH1")==NULL && strstr(s,"CH1CH1CH1")==NULL && k==13 && strchr(s,'2')==NULL && strstr(s,"CC4")==NULL){
	if(strstr(s,"CN1")||strstr(s,"CO1")||strstr(s,"CP1")||strstr(s,"CS1")||strstr(s,"Cl1")||strstr(s,"CF1")||strstr(s,"CB1")||strstr(s,"CI1")){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC3  ");
	printf("C3\n");}}	
	
	if(strstr(s,"CH1") && strchr(s,'2')==NULL && strstr(s,"CH1CH1CH1CH1")==NULL && strstr(s,"CH1CH1CH1")==NULL && strstr(s,"CH1CH1")==NULL && k==13 && strstr(s,"CC4")==NULL){
	if(strstr(s,"CN1")||strstr(s,"CO1")||strstr(s,"CP1")||strstr(s,"CS1")||strstr(s,"Cl1")||strstr(s,"CF1")||strstr(s,"CB1")||strstr(s,"CI1")){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC4  ");
	printf("C4\n");}}
	if(strstr(s,"CH1")==NULL && strchr(s,'2')==NULL && k==13 && strchr(s,'4')==NULL){
	if(strstr(s,"CN1")||strstr(s,"CO1")||strstr(s,"CP1")||strstr(s,"CS1")||strstr(s,"Cl1")||strstr(s,"CF1")||strstr(s,"CB1")||strstr(s,"CI1")){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC4  ");
	printf("C4\n");}}	
	
	if(strstr(s,"CO2") || strstr(s,"CN2") || strstr(s,"CS2") || strstr(s,"CP2")){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC5  ");
	printf("C5\n");}

	if(strstr(s,"CC2") && strstr(s,"CH1CH1") && k==10){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC6  ");
	printf("C6\n");}
	if(strstr(s,"CC2") && strstr(s,"CH1") && strstr(s,"CH1CH1")==NULL && strstr(s,"CC4")==NULL && k==10){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC6  ");
	printf("C6\n");}	
	if(strstr(s,"CC2") && strstr(s,"CH1")==NULL && strstr(s,"CC4")==NULL && k==10){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC6  ");
	printf("C6\n");}
	if(strstr(s,"CC2CC2")){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC6  ");
	printf("C6\n");}
	
        if(strchr(s,'3') && k==7){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC7  ");
	printf("C7\n");}

	if(strstr(s,"CH1CH1CH1") && strstr(s,"CC1") && k==16){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\tC8  ");
	printf("C8\n");}

	if(strstr(s,"CH1CH1CH1") && strstr(s,"CC1")==NULL && k==16){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\tC9  ");
	printf("C9\n");}

	if(strstr(s,"CH1CH1") && strstr(s,"CH1CH1CH1")==NULL && k==16){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\tC10  ");
	printf("C10\n");}

	if(strstr(s,"CH1") && strstr(s,"CH1CH1CH1")==NULL && strstr(s,"CH1CH1")==NULL && k==16){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\tC11  ");
	printf("C11\n");}

	if(strstr(s,"CH1")==NULL && strstr(s,"CH1CH1CH1")==NULL && strstr(s,"CH1CH1")==NULL && k==16){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\tC12  ");
	printf("C12\n");}

	if(strchr(s,'4') && strstr(s,"CP1") && k==10){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC13  ");
	printf("C13\n");}

	if(strchr(s,'4') && strstr(s,"CF1") && k==10){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC14  ");
	printf("C14\n");}
	if(strchr(s,'4') && strstr(s,"Cl1") && k==10){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC15  ");
	printf("C15\n");}
	if(strchr(s,'4') && strstr(s,"CB1") && k==10){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC16  ");
	printf("C16\n");}
	if(strchr(s,'4') && strstr(s,"CI1") && k==10){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC17  ");
	printf("C17\n");}

	if(strchr(s,'4') && strstr(s,"CH1") && k==10){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC18  ");
	printf("C18\n");}
	
	if(strchr(s,'4') && strchr(s,'1')==NULL && strchr(s,'2')==NULL && strchr(s,'3')==NULL && k==10){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC19  ");
	printf("C19\n");}

	if(strchr(s,'4') && strstr(s,"CC1") && strstr(s,"CX4") && k==13){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC20  ");
	printf("C20\n");}
		
	if(strchr(s,'4') && strstr(s,"CC1") && strstr(s,"CMe") && k==13){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC21  ");
	printf("C21\n");}

	if(strchr(s,'4') && strstr(s,"CN1") && k==10){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC22  ");
	printf("C22\n");}

	if(strchr(s,'4') && strstr(s,"CO1") && k==10){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC23  ");
	printf("C23\n");}

	if(strchr(s,'4') && strstr(s,"CS1") && k==10){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC24  ");
	printf("C24\n");}

	if(strstr(s,"CC2") && strstr(s,"CC4") && k==13){	
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC26  ");
	printf("C26\n");}

	if(strstr(s,"CC4CC4") || strstr(s,"CC4CN4") || strstr(s,"CC4CO4") || strstr(s,"CC4CS4") || strstr(s,"CC4CP4") || strstr(s,"CN4CP4") || strstr(s,"CN4CS4") || strstr(s,"CN4CO4") || strstr(s,"CN4CN4") || strstr(s,"CO4CO4") || strstr(s,"CO4CS4") || strstr(s,"CO4CP4") || strstr(s,"CS4CS4") || strstr(s,"CS4CP4") || strstr(s,"CP4CP4") || strstr(s,"CC4CC4")){
	if(k==10){
	if(strstr(s,"CO2") || strstr(s,"CC2") || strstr(s,"CN2"))
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tC25  ");
	printf("C25\n");}}	
	}
 	/********************************************************/	
	if(strncmp(s,"H",1)==0)
	{
	if(strstr(s,"HC1")){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\t\tH1  ");
	printf("H1\n");}

	if(strstr(s,"OOL") || strstr(s,"HP1") || strstr(s,"HS1") || strstr(s,"HF1") || strstr(s,"Hl1") || strstr(s,"HB1") || strstr(s,"HI1")){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\t\tH2  ");
	printf("H2\n");}
	
	if(strstr(s,"HN1") || strstr(s,"AMN")){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\t\tH3  ");
	printf("H3\n");}

	if(strstr(s,"ACD")){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\t\tH4  ");
	printf("H4\n");}
	}
	
	/*****************************************************/
	if(strncmp(s,"F",1)==0){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\t\tF  ");
	printf("F\n");}
	if(strncmp(s,"l",1)==0){
        strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\t\tCl  ");
	printf("Cl\n");}
	if(strncmp(s,"B",1)==0){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\t\tBr  ");
	printf("Br\n");}
	if(strncmp(s,"I",1)==0){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\t\tI  ");
	printf("I\n");}	
	
	/*********************************************************/
	//if(strncmp(s,"P",1)==0){
	//strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tP");
	//printf("P\n");}
	/*******New addition**********/
	if(strncmp(s,"P",1)==0){
	if(strstr(s,"PO2PO2PO2PO2")&&strstr(s,"OP2OP2OP2OP2"))
	strcat(s,"<--\t\t\t\t\t\t\t\t\tOP7");
        printf("OP7\n");
	if(strstr(s,"PO2PO2PO2")&&strstr(s,"PO2PO2PO2PO2")==NULL&&strstr(s,"OP2OP2OP2")&&strstr(s,"OP2OP2OP2OP2")==NULL)
	strcat(s,"<--\t\t\t\t\t\t\t\t\tOP8");
        printf("OP8\n");
	if(strstr(s,"PO2PO2PO1PO2")&&strstr(s,"PO2PO2PO2PO2")==NULL&&strstr(s,"OP2OP2OP1OP2")&&strstr(s,"OP2OP2OP2OP2")==NULL &&strstr(s,"OP2OP1OP2OP2")==NULL)
	strcat(s,"<--\t\t\t\t\t\t\t\t\tOP8");
        printf("OP8\n");
	if(strstr(s,"PO2PO1PO2PO2")&&strstr(s,"PO2PO2PO2PO2")==NULL&&strstr(s,"OP2OP1OP2OP2")&&strstr(s,"OP2OP2OP2OP2")==NULL && strstr(s,"OP2OP2OP1OP2")==NULL)
	strcat(s,"<--\t\t\t\t\t\t\t\t\tOP8");
        printf("OP8\n");
	if((strstr(s,"PO2PO2")||strstr(s,"PO2PO1PO2"))&&(strstr(s,"OP2OP2")||strstr(s,"OP2OP1OP2"))&&strstr(s,"PO2PO1PO2PO2")==NULL&&strstr(s,"PO2PO2PO1PO2")==NULL&&strstr(s,"PO2PO2PO2")==NULL&&strstr(s,"PO2PO2PO2PO2")==NULL&&strstr(s,"OP2OP2OP2OP2")==NULL&&strstr(s,"OP2OP1OP2OP2")==NULL&&strstr(s,"OP2OP2OP1OP2")==NULL&&strstr(s,"OP2OP2OP2")==NULL)
	strcat(s,"<--\t\t\t\t\t\t\t\t\tOP9");
        printf("OP9\n");}

	/*********************************************************/
	if(strncmp(s,"S",1)==0)
	{	
	if(k==7 && strchr(s,'2')==NULL && strchr(s,'4')==NULL){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tS1  ");
	printf("S1\n");}
	if(k==4 && strchr(s,'2')){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\t\tS1  ");
        printf("S1\n");}
	if(strchr(s,'2') && strchr(s,'4')==NULL && (k==10 || k==13)){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tS1  ");
	printf("S1\n");}
	if(strchr(s,'4')){
	strcat(s,"<--\t\t\t\t\t\t\t\t\t\t\tS3  ");
	printf("S3\n");}
	}
	fprintf(fw,"%s\n",s);
	}		
        printf("%s\n",s);	
	
	

     
		fclose(fr);
		fclose(fw);

}



