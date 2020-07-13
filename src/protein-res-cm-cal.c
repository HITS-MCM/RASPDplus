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
 *Purpose: Convert PDB format of a protein to RASPD(+) compatible format.
 *Please note that protein should not have any chain information and hydrogen added.
 * How to run:
 * ./protein-res-cm-cal.exe <PDBID of hydrogen added protein coordinates> 
 * say,
 * ./protein-res-cm-cal.exe 1NHZ-hydrogen.pdb
 *Input protein (protein.pdb):
ATOM      1  N   PRO     1       1.416  42.867 -14.065  1.00  0.00
ATOM      2  H2  PRO     1       2.203  42.394 -14.486  1.00  0.00
ATOM      3  H3  PRO     1       1.478  43.850 -13.843  1.00  0.00
ATOM      4  CD  PRO     1       0.236  42.766 -14.946  1.00  0.00
ATOM      5  HD2 PRO     1       0.572  42.609 -15.971  1.00  0.00
ATOM      6  HD3 PRO     1      -0.330  43.696 -14.888  1.00  0.00
ATOM      7  CG  PRO     1      -0.554  41.559 -14.389  1.00  0.00
ATOM      8  HG2 PRO     1      -0.712  40.842 -15.194  1.00  0.00
ATOM      9  HG3 PRO     1      -1.518  41.910 -14.021  1.00  0.00
ATOM     10  CB  PRO     1       0.262  40.961 -13.272  1.00  0.00
ATOM     11  HB2 PRO     1       0.844  40.120 -13.648  1.00  0.00
ATOM     12  HB3 PRO     1      -0.357  40.630 -12.438  1.00  0.00
ATOM     13  CA  PRO     1       1.203  42.082 -12.813  1.00  0.00
ATOM     14  HA  PRO     1       0.740  42.749 -12.086  1.00  0.00
ATOM     15  C   PRO     1       2.482  41.507 -12.229  1.00  0.00
ATOM     16  O   PRO     1       3.293  40.917 -12.935  1.00  0.00
ATOM     17  N   THR     2       2.666  41.692 -10.930  1.00  0.00
ATOM     18  H   THR     2       2.017  42.249 -10.393  1.00  0.00
ATOM     19  CA  THR     2       3.756  41.027 -10.224  1.00  0.00
ATOM     20  HA  THR     2       4.685  41.355 -10.691  1.00  0.00
ATOM     21  CB  THR     2       3.788  41.412  -8.762  1.00  0.00
ATOM     22  HB  THR     2       4.665  40.954  -8.305  1.00  0.00
ATOM     23  CG2 THR     2       3.853  42.961  -8.575  1.00  0.00
ATOM     24 HG21 THR     2       2.976  43.420  -9.031  1.00  0.00
ATOM     25 HG22 THR     2       3.875  43.199  -7.511  1.00  0.00
ATOM     26 HG23 THR     2       4.754  43.348  -9.051  1.00  0.00
ATOM     27  OG1 THR     2       2.556  40.968  -8.164  1.00  0.00
ATOM     28  HG1 THR     2       2.550  41.201  -7.233  1.00  0.00
ATOM     29  C   THR     2       3.549  39.531 -10.221  1.00  0.00
ATOM     30  O   THR     2       2.427  39.020 -10.347  1.00  0.00
*Output from this program (./protein-res-cm-cal.exe protein.pdb):
PRO     O       1            3.293   40.917  -12.935
PRO     CM      1            1.202   41.819  -13.558
THR     H       2            2.017   42.249  -10.393
THR     OG1     2            2.556   40.968   -8.164
THR     HG1     2            2.550   41.201   -7.233
THR     O       2            2.427   39.020  -10.347
THR     CM      2            3.187   40.974   -9.552 
*/

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
int main(int argc, char *argv[])
{
        FILE *fr;
	if (argc<2) {
  	printf("./protein-res-cm-cal.exe <PDBID of hydrogen added protein coordinates>\n");
  	exit(1);
	}
        char str[4], s[100];
	float sum=0, m=0;
        float d1=0,d2=0,d3=0;
        float mx=0,my=0,mz=0;
        struct format
        {
                char atom[6], atmsmb[6], resid[6];
                int atomno, resno;
                float x, y, z, chg, rad,epsi;
        };
        struct format f;

        fr=fopen(argv[1],"r");
        int start=1,last;
        while(fgets(s,100,fr)!=NULL)
        {
                sscanf(s,"%s%d%s%s%d%f%f%f%f%f%f",f.atom,&f.atomno,f.atmsmb,f.resid,&f.resno,&f.x,&f.y,&f.z,&f.chg,&f.rad,&f.epsi);
		if((strchr(f.atmsmb,'O') && strlen(f.atmsmb)==1) || (strchr(f.atmsmb,'H') && strlen(f.atmsmb)==1) ||
		(strstr(s,"ARG") && strstr(f.atmsmb,"1HH1")) ||
		(strstr(s,"ARG") && strstr(f.atmsmb,"HH11")) ||
		(strstr(s,"ARG") && strstr(f.atmsmb,"2HH1")) ||
		(strstr(s,"ARG") && strstr(f.atmsmb,"HH12")) ||
		(strstr(s,"ARG") && strstr(f.atmsmb,"1HH2")) ||
		(strstr(s,"ARG") && strstr(f.atmsmb,"HH21")) ||
		(strstr(s,"ARG") && strstr(f.atmsmb,"2HH2")) ||
		(strstr(s,"ARG") && strstr(f.atmsmb,"HH22")) ||
		(strstr(s,"ARG") && strstr(f.atmsmb,"HE"))   ||
		(strstr(s,"LYS") && strstr(f.atmsmb,"HZ"))   ||
		(strstr(s,"LYN") && strstr(f.atmsmb,"HZ"))   ||
		(strstr(s,"ASN") && strstr(f.atmsmb,"1HD2")) ||
		(strstr(s,"ASN") && strstr(f.atmsmb,"HD21")) ||
		(strstr(s,"ASN") && strstr(f.atmsmb,"2HD2")) ||
                (strstr(s,"ASN") && strstr(f.atmsmb,"HD22")) ||
		(strstr(s,"GLN") && strstr(f.atmsmb,"1HE2")) ||
		(strstr(s,"GLN") && strstr(f.atmsmb,"HE21")) ||
		(strstr(s,"GLN") && strstr(f.atmsmb,"2HE2")) ||
		(strstr(s,"GLN") && strstr(f.atmsmb,"HE22")) ||
		(strstr(s,"THR") && strstr(f.atmsmb,"1HG"))  ||
		(strstr(s,"THR") && strstr(f.atmsmb,"HG1"))  ||
		(strstr(s,"SER") && strstr(f.atmsmb,"HG"))   ||
		(strstr(s,"TYR") && strstr(f.atmsmb,"HH"))   ||
		(strstr(s,"HIE") && strstr(f.atmsmb,"HE2"))  ||
		(strstr(s,"HID") && strstr(f.atmsmb,"HD1"))  ||
		(strstr(s,"HIP") && strstr(f.atmsmb,"HD1"))  ||
		(strstr(s,"HIP") && strstr(f.atmsmb,"HE2"))  ||
		(strstr(s,"ASH") && strstr(f.atmsmb,"HD"))   ||
		(strstr(s,"GLH") && strstr(f.atmsmb,"HE"))   ||
		(strstr(s,"TRP") && strstr(f.atmsmb,"HE1"))  ||
		(strstr(s,"TRP") && strstr(f.atmsmb,"1HE"))  ||
		(strstr(s,"THR") && strstr(f.atmsmb,"OG1"))  ||
		(strstr(s,"SER") && strstr(f.atmsmb,"OG"))   ||
		(strstr(s,"TYR") && strstr(f.atmsmb,"OH"))   ||
		(strstr(s,"ASN") && strstr(f.atmsmb,"OD1"))  ||
		(strstr(s,"GLN") && strstr(f.atmsmb,"OE1"))  ||
		(strstr(s,"ASP") && strstr(f.atmsmb,"OD"))   ||
		(strstr(s,"GLU") && strstr(f.atmsmb,"OE"))   ||
		(strstr(s,"HIE") && strstr(f.atmsmb,"ND"))   ||
		(strstr(s,"LYN") && strstr(f.atmsmb,"NZ"))   ||
		(strstr(s,"HID") && strstr(f.atmsmb,"NE"))   ||
		(strstr(s,"ASH") && strstr(f.atmsmb,"OD"))   ||
		(strstr(s,"GLH") && strstr(f.atmsmb,"OE")))

		printf("%s\t%s\t%d\t  %8.3f %8.3f %8.3f\n", f.resid, f.atmsmb, f.resno, f.x, f.y, f.z);
		if (f.resno==start){
			strcpy(str,f.resid);
				if(strncmp(f.atmsmb,"C",1)==0 || strncmp(f.atmsmb,"1C",2)==0 || strncmp(f.atmsmb,"2C",2)==0 || strncmp(f.atmsmb,"3C",2)==0)
                                        m=12;
                                if(strncmp(f.atmsmb,"N",1)==0 || strncmp(f.atmsmb,"1N",2)==0 || strncmp(f.atmsmb,"2N",2)==0 || strncmp(f.atmsmb,"3N",2)==0)
                                        m=14;
                                if(strncmp(f.atmsmb,"H",1)==0 || strncmp(f.atmsmb,"1H",2)==0 || strncmp(f.atmsmb,"2H",2)==0 || strncmp(f.atmsmb,"3H",2)==0)
                                        m=1;
                                if(strncmp(f.atmsmb,"O",1)==0 || strncmp(f.atmsmb,"1O",2)==0 || strncmp(f.atmsmb,"2O",2)==0 || strncmp(f.atmsmb,"3O",2)==0)
                                        m=16;
                                if(strncmp(f.atmsmb,"S",1)==0 || strncmp(f.atmsmb,"1S",2)==0 || strncmp(f.atmsmb,"2S",2)==0 || strncmp(f.atmsmb,"3S",2)==0)
                                        m=32;
                                if(strncmp(f.atmsmb,"Se",2)==0 || strncmp(f.atmsmb,"1Se",3)==0 || strncmp(f.atmsmb,"2Se",3)==0 || strncmp(f.atmsmb,"3Se",3)==0)
                                        m=79;
                                if(strncmp(f.atmsmb,"P",1)==0 || strncmp(f.atmsmb,"1P",2)==0 || strncmp(f.atmsmb,"2P",2)==0 || strncmp(f.atmsmb,"3P",2)==0)
                                        m=31;
				sum=sum+m;
                                mx=mx+m*f.x;
                                my=my+m*f.y;
                                mz=mz+m*f.z;
				last=f.resno;
                }
else {
        d1=mx/sum;
        d2=my/sum;
        d3=mz/sum;
        printf("%s\tCM\t%d\t  %8.3f %8.3f %8.3f\n", str, f.resno-1, d1, d2, d3);
        mx=my=mz=sum=0;
                                if(strncmp(f.atmsmb,"C",1)==0 || strncmp(f.atmsmb,"1C",2)==0 || strncmp(f.atmsmb,"2C",2)==0 || strncmp(f.atmsmb,"3C",2)==0)
                                        m=12;
                                if(strncmp(f.atmsmb,"N",1)==0 || strncmp(f.atmsmb,"1N",2)==0 || strncmp(f.atmsmb,"2N",2)==0 || strncmp(f.atmsmb,"3N",2)==0)
                                        m=14;
                                if(strncmp(f.atmsmb,"H",1)==0 || strncmp(f.atmsmb,"1H",2)==0 || strncmp(f.atmsmb,"2H",2)==0 || strncmp(f.atmsmb,"3H",2)==0)
                                        m=1;
                                if(strncmp(f.atmsmb,"O",1)==0 || strncmp(f.atmsmb,"1O",2)==0 || strncmp(f.atmsmb,"2O",2)==0 || strncmp(f.atmsmb,"3O",2)==0)
                                        m=16;
                                if(strncmp(f.atmsmb,"S",1)==0 || strncmp(f.atmsmb,"1S",2)==0 || strncmp(f.atmsmb,"2S",2)==0 || strncmp(f.atmsmb,"3S",2)==0)
                                        m=32;
                                if(strncmp(f.atmsmb,"Se",2)==0 || strncmp(f.atmsmb,"1Se",3)==0 || strncmp(f.atmsmb,"2Se",3)==0 || strncmp(f.atmsmb,"3Se",3)==0)
                                        m=79;
                                if(strncmp(f.atmsmb,"P",1)==0 || strncmp(f.atmsmb,"1P",2)==0 || strncmp(f.atmsmb,"2P",2)==0 || strncmp(f.atmsmb,"3P",2)==0)
                                        m=31;                                
				sum=sum+m;
                                mx=mx+m*f.x;
                                my=my+m*f.y;
                                mz=mz+m*f.z;

                start=f.resno;
                }
        }
        d1=mx/sum;
        d2=my/sum;
        d3=mz/sum;
        printf("%s\tCM\t%d\t  %8.3f %8.3f %8.3f\n", str, last, d1, d2, d3);
        fclose(fr);
}
