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
 * Purpose: This program fetch the values of logP and MR for the atom types of atoms in a molecule from the work of Scott A. Wildman and Gordon M. Crippen (J. Chem. Inf. Comput. Sci. 1999, 39, 868).
 * How to run:
 * ./lp8.exe <Output file of lp7.exe> <Output file name1> <Output file name2> 
*/

#include<stdio.h>
#include<string.h>
#include<stdlib.h>

int main(int argc, char *argv[])
{
        FILE *fr,*fw1,*fw2;
	if (argc<2) {
        printf("               <Input file-1>        <Output file-1>      <Output file-2>\n");
        printf("./lp8.exe <Output file of lp7.exe> <Output file name1> <Output file name2>\n");
        exit(1);
        }

        char strings[500],str[50],s[200];
        char strings_array[500][500];
        char str_array[500][500];
        char final_string[500][500];
        char *result;
        int k;
	float sum1=0,sum2=0;
        fr=fopen(argv[1],"r");
        fw1=fopen(argv[2],"w");
	fw2=fopen(argv[3],"w");

	if (fw1==NULL || fw2==NULL)            
        {
        printf("lp8.exe: Cannot write output. No output file name is given\n");
        exit(1);
        }

        int line=0;
        while(fgets(strings,500,fr)!=NULL)
        {
	//sum1=0;
	s[0]='\0';
        sscanf(strings,"%s%s",str,s);
	k=strlen(strings);
	if(strstr(s,"C1") &&strstr(s,"C10")==NULL&& strstr(s,"C11")==NULL&&strstr(s,"C12")==NULL&&strstr(s,"C13")==NULL&&strstr(s,"C14")==NULL&&strstr(s,"C15")==NULL&&strstr(s,"C16")==NULL&&strstr(s,"C17")==NULL&&strstr(s,"C18")==NULL&&strstr(s,"C19")==NULL){
	sum1=sum1+0.1441;
	sum2=sum2+2.503;}
	if(strstr(s,"C2") && strstr(s,"C20")==NULL&&strstr(s,"C21")==NULL&&strstr(s,"C22")==NULL&&strstr(s,"C23")==NULL&&strstr(s,"C24")==NULL&&strstr(s,"C25")==NULL&&strstr(s,"C26")==NULL){
        sum1=sum1+0.0000;
        sum2=sum2+2.433;}
	if(strstr(s,"C3")){
        sum1=sum1+-0.2035;
        sum2=sum2+2.753;}
	if(strstr(s,"C4")){
        sum1=sum1+-0.2051;
        sum2=sum2+2.731;}
	if(strstr(s,"C5")){
        sum1=sum1+-0.2783;
        sum2=sum2+5.007;}
	if(strstr(s,"C6")){
        sum1=sum1+0.1551;
        sum2=sum2+3.513;}
	if(strstr(s,"C7")){
        sum1=sum1+0.00170;
        sum2=sum2+3.888;}
	if(strstr(s,"C8")){
        sum1=sum1+0.08452;
        sum2=sum2+2.464;}
	if(strstr(s,"C9")){
        sum1=sum1+-0.1444;
        sum2=sum2+2.412;}
	if(strstr(s,"C10")){
        sum1=sum1+-0.0516;
        sum2=sum2+2.488;}
	if(strstr(s,"C11")){
        sum1=sum1+0.1193;
        sum2=sum2+2.582;}
	if(strstr(s,"C12")){
        sum1=sum1+-0.0967;
        sum2=sum2+2.576;}
	if(strstr(s,"C13")){
        sum1=sum1+-0.5443;
        sum2=sum2+4.041;}
	if(strstr(s,"C14")){
        sum1=sum1+0.0000;
        sum2=sum2+3.257;}
	if(strstr(s,"C15")){
        sum1=sum1+0.2450;
        sum2=sum2+3.564;}
	if(strstr(s,"C16")){
        sum1=sum1+0.1980;
        sum2=sum2+3.180;}
	if(strstr(s,"C17")){
        sum1=sum1+0.0000;
        sum2=sum2+3.104;}
	if(strstr(s,"C18")){
        sum1=sum1+0.1581;
        sum2=sum2+3.350;}
	if(strstr(s,"C19")){
        sum1=sum1+0.2955;
        sum2=sum2+4.346;}
	if(strstr(s,"C20")){
        sum1=sum1+0.2713;
        sum2=sum2+3.904;}
	if(strstr(s,"C21")){
        sum1=sum1+0.1360;
        sum2=sum2+3.509;}
	if(strstr(s,"C22")){
        sum1=sum1+0.4619;
        sum2=sum2+4.067;}
	if(strstr(s,"C23")){
        sum1=sum1+0.5437;
        sum2=sum2+3.853;}
	if(strstr(s,"C24")){
        sum1=sum1+0.1893;
        sum2=sum2+2.673;}
	if(strstr(s,"C25")){
        sum1=sum1+-0.8186;
        sum2=sum2+3.135;}
	if(strstr(s,"C26")){
        sum1=sum1+0.2640;
        sum2=sum2+4.305;}
	if(strstr(str,"<-")==NULL){
	if(strncmp(str,"C",1)==0){
	sum1=sum1+0.08129;
	sum2=sum2+3.243;}}
	
	if(strstr(s,"N1") && strstr(s,"N10")==NULL&& strstr(s,"N11")==NULL&&strstr(s,"N12")==NULL&&strstr(s,"N13")==NULL&&strstr(s,"N14")==NULL&&strstr(s,"N15")==NULL&&strstr(s,"N16")==NULL&&strstr(s,"N17")==NULL&&strstr(s,"N18")==NULL&&strstr(s,"N19")==NULL){
	sum1=sum1+-1.0190;
	sum2=sum2+2.262;}
	if(strstr(s,"N2")){
        sum1=sum1+-0.7096;
        sum2=sum2+2.173;}
	if(strstr(s,"N3")){
        sum1=sum1+-1.0270;
        sum2=sum2+2.827;}
	if(strstr(s,"N4")){
        sum1=sum1+-0.5188;
        sum2=sum2+3.000;}
	if(strstr(s,"N5")){
        sum1=sum1+0.08387;
        sum2=sum2+1.757;}
	if(strstr(s,"N6")){
        sum1=sum1+0.1836;
        sum2=sum2+2.428;}
	if(strstr(s,"N7")){
        sum1=sum1+-0.3187;
        sum2=sum2+1.839;}
	if(strstr(s,"N8")){
        sum1=sum1+-0.4458;
        sum2=sum2+2.819;}
	if(strstr(s,"N9")){
        sum1=sum1+0.01508;
        sum2=sum2+1.725;}
	if(strstr(s,"N10")){
        sum1=sum1+-1.950;
        sum2=sum2+0.0000;}
	if(strstr(s,"N11")){
        sum1=sum1+-0.3239;
        sum2=sum2+2.202;}
	if(strstr(s,"N12")){
        sum1=sum1+-1.119;
        sum2=sum2+0.0000;}
	if(strstr(s,"N13")){
        sum1=sum1+-0.3396;
        sum2=sum2+0.2604;}
	if(strstr(s,"N14")){
        sum1=sum1+0.2887;
        sum2=sum2+3.359;}
	if(strstr(str,"<-")==NULL){
        if(strncmp(str,"N",1)==0){
        sum1=sum1+-0.4806;
        sum2=sum2+2.134;}}

	if(strstr(s,"O1") &&strstr(s,"O10")==NULL&& strstr(s,"O11")==NULL&&strstr(s,"O12")==NULL&&strstr(s,"O13")==NULL&&strstr(s,"O14")==NULL&&strstr(s,"O15")==NULL&&strstr(s,"O16")==NULL&&strstr(s,"O17")==NULL&&strstr(s,"O18")==NULL&&strstr(s,"O19")==NULL){
	sum1=sum1+0.1552;
	sum2=sum2+1.080;}
	if(strstr(s,"O2")){
        sum1=sum1+-0.2893;
        sum2=sum2+0.8238;}
	if(strstr(s,"O3")){
        sum1=sum1+-0.0684;
        sum2=sum2+1.085;}
	if(strstr(s,"O4")){
        sum1=sum1+-0.4195;
        sum2=sum2+1.182;}
	if(strstr(s,"O5")){
        sum1=sum1+0.0335;
        sum2=sum2+3.367;}
	if(strstr(s,"O6")){
        sum1=sum1+-0.3339;
        sum2=sum2+0.7774;}
	if(strstr(s,"O7")){
        sum1=sum1+-1.189;
        sum2=sum2+0.000;}
	if(strstr(s,"O8")){
        sum1=sum1+0.1788;
        sum2=sum2+3.135;}
	if(strstr(s,"O9")){
        sum1=sum1+-0.1526;
        sum2=sum2+0.000;}
	if(strstr(s,"O10")){
        sum1=sum1+0.1129;
        sum2=sum2+0.2215;}
	if(strstr(s,"O11")){
        sum1=sum1+0.4833;
        sum2=sum2+0.3890;}
	if(strstr(s,"O12")){
        sum1=sum1+-1.326;
        sum2=sum2+0.0000;}
	if(strstr(str,"<-")==NULL){
        if(strncmp(str,"O",1)==0){
        sum1=sum1+-0.1188;
        sum2=sum2+0.6865;}}

	if(strchr(s,'F')){
	sum1=sum1+0.4202;
	sum2=sum2+1.108;}
	if(strstr(s,"Cl")){
        sum1=sum1+0.6895;
        sum2=sum2+5.853;}
	if(strstr(s,"Br")){
        sum1=sum1+0.8456;
        sum2=sum2+8.927;}
	if(strchr(s,'I')){
        sum1=sum1+0.8857;
        sum2=sum2+14.02;}
	
	if(strstr(s,"H1")){
        sum1=sum1+0.1230;
        sum2=sum2+1.057;}
	if(strstr(s,"H2")){
        sum1=sum1+-0.2677;
        sum2=sum2+1.395;}
	if(strstr(s,"H3")){
        sum1=sum1+0.2142;
        sum2=sum2+0.9627;}
	if(strstr(s,"H4")){
        sum1=sum1+0.2980;
        sum2=sum2+1.805;}
	if(strstr(str,"<-")==NULL){
        if(strncmp(str,"H",1)==0){
        sum1=sum1+0.1125;
        sum2=sum2+1.112;}}

	//if(strchr(s,'P')){
	//sum1=sum1+0.8612;
	//sum2=sum2+6.920;}
	if(strstr(s,"OP7")){
	sum1=sum1+-2.3494;  
	sum2=sum2+6.920;}
	if(strstr(s,"OP8")){
	sum1=sum1+-1.2792; 
        sum2=sum2+6.920;}
	if(strstr(s,"OP9")){
        sum1=sum1+-0.209;  
        sum2=sum2+6.920;}
	if(strstr(str,"<-")==NULL){
        if(strncmp(str,"P",1)==0){
	sum1=sum1+0.8612;
        sum2=sum2+6.920;}}

	if(strstr(s,"S1")){
	sum1=sum1+0.6482;
	sum2=sum2+7.591;}
	if(strstr(s,"S3")){
        sum1=sum1+0.6237;
        sum2=sum2+6.691;}
	if(strstr(str,"<-")==NULL){
        if(strncmp(str,"S",1)==0){
        sum1=sum1+-0.0024;
        sum2=sum2+7.365;}}

	//printf("%f\n",sum1);
}
	fprintf(fw1,"%f\n",sum1);
	fprintf(fw2,"%f\n",sum2);
fclose(fr);
}
