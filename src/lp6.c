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
 * Purpose: This program will concatenate Phosphorus (P) atom types with all the neighbouring group atom types which are attach with P atom
 * Say, in Fosinopril molecule, there is one P atom which is connected with two carbons atoms and one double bonded oxygen and one single bonded oxygen. Therefore, the primary atom types of P atom:
 * PC1PC1PO1PO2
 * This program will concatenate all the neighbouring group atom types which are attach with P atom. It means,
 * The atom types of two carbon atoms and two oxygen atoms will be concatenated with the central P atom by this program. The final atom type of P atom:
 * PC1PC1PO1PO2 CC1CP1CH1CH1 CC1CP1CH1CH1 OC1OP1 OP2
 * How to run:
 * ./lp6.exe <Output file of lp5.exe> <Output file of lp3.exe> <Output file name> 
*/

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
//#include <iostream.h>

int main(int argc, char *argv[])
{
	FILE *fp,*fr,*fw;
	if (argc<3) {
  	printf("               <Input file-1>            <Input file-2>       <Output file>\n");
  	printf("./lp6.exe <Output file of lp5.exe> <Output file of lp3.exe> <Output file name>\n");
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

        char strings[100],str[100];
	char strings_array[1000][100];
	//char str_array[100][100];
	char final_string[300][100];
 	char *result;	
	fp=fopen(argv[1],"r");
	fr=fopen(argv[2],"r");
	fw=fopen(argv[3],"w");
		
	if (fw==NULL)            
        {
        printf("lp6.exe: Cannot write output. No output file name is given\n");
        exit(1);
        }

	int line1=0;
	while(fgets(strings,100,fp)!=NULL)
	{
		int len = strlen(strings);
		strncpy(strings_array[line1],strings,len-1);	
		line1++;
	}


	int line=0;
	while(fgets(str,100,fr)!=NULL)
	{
                f.a=0;
                f.b=0;
                f.c=0;
		f.d=0;
		f.e=0;
		f.f=0;
		f.g=0;
                sscanf(str,"%d%d%d%d%d%d%d",&f.a,&f.b,&f.c,&f.d,&f.e,&f.f,&f.g);
                    
                
		strcpy(final_string[line],strings_array[line]);
                if(strncmp(final_string[line],"P",1)==0){       	                
	        strncat(final_string[line],strings_array[f.b-1],18);          
	        strncat(final_string[line],strings_array[f.c-1],18);
                strncat(final_string[line],strings_array[f.d-1],18);
		strncat(final_string[line],strings_array[f.e-1],18);
		strncat(final_string[line],strings_array[f.f-1],18);}
		//strncat(final_string[line],strings_array[f.g-1],18);
                fprintf(fw,"%s\n",final_string[line]);
		
		

		line++;
	        
              
              //printf("%s\n",final_string[line]);
          }
          

	fclose(fp);
	fclose(fr);

}

