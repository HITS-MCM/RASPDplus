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

 * Written by Satyanarayana Rao (SCFBio, IITD currently at University of Colorado, Denve) and Goutam Mukherjee
 *Purpose: Calculate physicochemical parameters of active site of a protein.
 * How to run:
 * ./protein_param_calc <Output of protein-res-cm-cal.exe> <Output of maxD_calculator.exe> <Center of mass of protein active site (this information is passed as an argument)> 
 * Say, protein.ent file (Output of protein-res-cm-cal.exe) contains the following input:
 * ASP     H       109         -2.313    5.369   17.030
 * ASP     OD1     109         -1.498    3.169   14.916
 * ASP     OD2     109         -2.434    2.870   12.991
 * ASP     O       109         -1.283    7.866   13.749
 * ASP     CM      109         -1.988    5.160   14.364
 * MET     H       110         -3.464    7.840   15.995
 * MET     O       110         -1.338   11.282   14.555
 * MET     CM      110         -4.043    9.982   14.225
 * LYS     H       115          2.005   13.450   14.681
 * LYS     HZ1     115          6.441   10.205   17.653
 * LYS     HZ2     115          7.777   10.143   16.687
 * LYS     HZ3     115          7.595   11.378   17.765
 * LYS     O       115          4.709   16.805   14.163
 * LYS     CM      115          4.677   13.804   15.502
 * HIE     H       116          3.786   14.156   12.736
 * HIE     ND1     116          7.280   13.139   12.013
 * HIE     HE2     116          6.448   10.173   12.448
 * HIE     O       116          5.305   17.581   10.562
 * HIE     CM      116          5.761   14.059   11.671
 * MET     H       117          2.735   15.871   11.084
 * MET     O       117          2.306   20.118   10.386
 * MET     CM      117          1.372   16.730    9.423
 *
 * Dmax value (Output of maxD_calculator.exe) is 8.738 Å (say)
 * Please not here, in Dmax file, values should be written twice (say, 8.738 8.738). Since, this program takes first value to set distance cut off and the next value for normalization. In RASPD/RASPD+, both the values are same, hence, we have to paste the Dmax values twice in the file.
 *
 * Center of mass of protein active site is -2.251  12.432   3.402
 * Then,
 * ./protein_param_calc protein.ent Dmax -2.251  12.432   3.402
 * Output of this program:
 * 0.000	 0.000	 0.000	 0.000	 0.000	 0.000	 0.000	 0.000	 0.229	 0.458	 0.201	 0.000	 4.419	 0.000
 * Column-1: PD(K+R+HIP); Column-2: PA(D+E); Column-3: PD(T+S+Y+DH+EH); Column-4: PA(N+Q+T+S+DH+EH); Column-5: PD(LYN+N+Q); Column-6: PA(LYN); Column-7: PD(W+H); Column-8: PA(Y+H); Column-9: PD(Amide-NH);Column-10: PA(Amide-O); Column-11: PP(Non-Arom); Column-12: PP(Arom);Column-13: PMR(Non-Arom); Column-14: PMR(Arom)
 
 * PA: Protein Acceptor
 * PD: Protein Donor
 * PP: logP values of protein active site residues
 * PMR: Molar refractivities of protein active site residues
 * One letter code of amino acid code is given in parentheses
 * DH: Hydrogen added Aspartic acid
 * EH: Hydrogen added Glutamic acid
 * LYN: Neutral Lysine residue (Never appear in training set)
 * Non-Arom: Non aromatic residues
 * Arom: Aromatic residues
*/

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include <math.h>
#include <cfloat>
#include <ctime>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <map>


using namespace std;
map <string, int> HBD,HBA;
map <string, float> logP1;
map <string, float> logP2;
map <string, float> pmr1;
map <string, float> pmr2;
map <string, int>::iterator it_HBD, it_HBA;
map <string, float>::iterator it_logP1;
map <string, float>::iterator it_pmr1;
map <string, float>::iterator it_logP2;
map <string, float>::iterator it_pmr2;

vector <float> maxdistance;
vector <float> rog_map;
void intiateHbond ();
typedef struct ligand {
	float x,y,z;
	int r;
	char  atomname[5];
	char resname[4];
}L;
typedef struct Protein {
	float x,y,z;
	int r;
	char  atomname[5];
	char  resname[4];

}P;

int main(int argc, char *argv[])
{
	FILE *fp, *fd;
	if (argc<6) {
	 printf("                 	 <Input file-1> <Input file-2>  <Input arguments>\n");	
 	 printf("./protein_param_calc.exe  protein.ent      Dmax      -2.251  12.432   3.402\n");
 	 exit(1);
	}
        char s[100];
	time_t start, end;
	float mdis=0,rog=0;
	
	time(&start);
	fp=fopen(argv[1],"r");
	fd=fopen(argv[2],"r");

	float com_x=atof(argv[3]);
	float com_y=atof(argv[4]);
	float com_z=atof(argv[5]);
	
	int counter_protein=0, counter_drug=0, donor=0, acceptor=0;
	int sum_HBA=0, sum_HBD=0;
	int sum_HBD1=0, sum_HBD2=0;
	int sum_HBD3=0;
	int sum_HBD4=0;
        int sum_HBA1=0, sum_HBA2=0;
        int sum_HBA3=0;
        int sum_HBA4=0;

        float sum_logP1=0, sum_pmr1=0, sum_logP2=0, sum_pmr2=0;

	P Pro [200000];
	L Lig [500];
	vector <string> PDB;
	intiateHbond();
	while (fgets (s, 100, fp)!=NULL){
		if (6==sscanf(s,"%s%s%d%f%f%f",Pro[counter_protein].atomname\
					,Pro[counter_protein].resname,&Pro[counter_protein].r,&Pro[counter_protein].x\
					,&Pro[counter_protein].y,&Pro[counter_protein].z)){
			PDB.push_back(s);
			counter_protein++;
		}
	}
	while (fgets (s, 100, fd)!=NULL){
		sscanf(s,"%f %f", &mdis, &rog);
		maxdistance.push_back(mdis);
		rog_map.push_back(rog);
	}
	double distance=0;
	for (int i=0; i<maxdistance.size(); i++){
		sum_HBD=sum_HBA=0;
                sum_HBD1=sum_HBD2=sum_HBD3=sum_HBD4=0;
                sum_HBA1=sum_HBA2=sum_HBA3=sum_HBA4=0;
		sum_logP1=sum_pmr1=0;
		sum_logP2=sum_pmr2=0;
		for (int j =0; j<counter_protein;j++){
			distance=sqrt((Pro[j].x - com_x)*(Pro[j].x - com_x) + (Pro[j].y - com_y)*(Pro[j].y - com_y) + (Pro[j].z - com_z)*(Pro[j].z - com_z));
			if (distance <= maxdistance.at(i) + 3){
				if (strcmp(Pro[j].resname, "CM")){
					//Charged HBD and HBA count at the protein active site
					string tmp(Pro[j].atomname);
					tmp+="-";
					tmp.append(Pro[j].resname);
					if (HBD.find(tmp)!=HBD.end()){
						sum_HBD++;
					}
					if (HBA.find(tmp)!=HBA.end()){
                                                sum_HBA++;
                                        }
					//Neutral N-containing HBD and HBA count at the protein active site
					string tmp1(Pro[j].atomname);
                                        tmp1+="+";
                                        tmp1.append(Pro[j].resname);
                                        if (HBD.find(tmp1)!=HBD.end()){
                                                sum_HBD1++;
                                        }
					if (HBA.find(tmp1)!=HBA.end()){
                                                sum_HBA1++;
                                        }
					//Neutral O-containing HBD and HBA count at the protein active site
					string tmp2(Pro[j].atomname);
                                        tmp2+="*";
                                        tmp2.append(Pro[j].resname);
                                        if (HBD.find(tmp2)!=HBD.end()){
                                                sum_HBD2++;
                                        }
					if (HBA.find(tmp2)!=HBA.end()){
                                                sum_HBA2++;
                                        }
					//Aromatic HBD and HBA count at the protein active site
					string tmp3(Pro[j].atomname);
                                        tmp3+="#";
                                        tmp3.append(Pro[j].resname);
                                        if (HBD.find(tmp3)!=HBD.end()){
                                                sum_HBD3++;
                                        }
					if (HBA.find(tmp3)!=HBA.end()){
                                                sum_HBA3++;
                                        }
					//Amide group of peptide bonds
					if (!strcmp(Pro[j].resname,"H")){
                                                sum_HBD4++;
                                        }
					if (!strcmp(Pro[j].resname,"O")){
                                                sum_HBA4++;
                                        }
				}
					else if (distance <= maxdistance.at(i) + 0.9){
					if (pmr1.find(Pro[j].atomname)!=pmr1.end()){
						sum_pmr1 +=pmr1[Pro[j].atomname];
					}
					if (pmr2.find(Pro[j].atomname)!=pmr2.end()){
                                                sum_pmr2 +=pmr2[Pro[j].atomname];
                                        }
					if (logP1.find(Pro[j].atomname)!=logP1.end()){
						sum_logP1 +=logP1[Pro[j].atomname];
					}
					if (logP2.find(Pro[j].atomname)!=logP2.end()){
                                                sum_logP2 +=logP2[Pro[j].atomname];
                                        }

				}
			}

		}
		//cout<<sum_HBD/rog_map.at(i) <<"\t"<<sum_HBA/rog_map.at(i) <<"\t"<<sum_logP/rog_map.at(i)<<"\t"<<sum_pmr/rog_map.at(i)<<"\t"<<"7.207"<<endl;
		printf("%0.3f\t %0.3f\t %0.3f\t %0.3f\t %0.3f\t %0.3f\t %0.3f\t %0.3f\t %0.3f\t %0.3f\t %0.3f\t %0.3f\t %0.3f\t %0.3f\n", sum_HBD/rog_map.at(i), sum_HBA/rog_map.at(i),  sum_HBD2/rog_map.at(i), sum_HBA2/rog_map.at(i), sum_HBD1/rog_map.at(i), sum_HBA1/rog_map.at(i), sum_HBD3/rog_map.at(i), sum_HBA3/rog_map.at(i), sum_HBD4/rog_map.at(i), sum_HBA4/rog_map.at(i), sum_logP1/rog_map.at(i), sum_logP2/rog_map.at(i), sum_pmr1/rog_map.at(i), sum_pmr2/rog_map.at(i));
		//printf("%0.3f\t %0.3f\t %0.3f\t %0.3f\n", sum_logP1/rog_map.at(i), sum_logP2/rog_map.at(i), sum_pmr1/rog_map.at(i), sum_pmr2/rog_map.at(i));
	}
	

}
void intiateHbond (){
//Charged-HBD
	HBD.insert(pair<string,int>("ARG-1HH1",1));
        HBD.insert(pair<string,int>("ARG-2HH1",1));
        HBD.insert(pair<string,int>("ARG-1HH2",1));
        HBD.insert(pair<string,int>("ARG-2HH2",1));
        HBD.insert(pair<string,int>("ARG-HE",1));  //Additional line as the format of tleap
        HBD.insert(pair<string,int>("ARG-HH11",1)); //Additional line as the format of tleap
        HBD.insert(pair<string,int>("ARG-HH12",1)); //Additional line as the format of tleap
        HBD.insert(pair<string,int>("ARG-HH21",1)); //Additional line as the format of tleap
        HBD.insert(pair<string,int>("ARG-HH22",1)); //Additional line as the format of tleap
        HBD.insert(pair<string,int>("LYS-1HZ",1));
        HBD.insert(pair<string,int>("LYS-2HZ",1));
        HBD.insert(pair<string,int>("LYS-3HZ",1));
        HBD.insert(pair<string,int>("LYS-HZ1",1)); //Additional line as the format of tleap
        HBD.insert(pair<string,int>("LYS-HZ2",1)); //Additional line as the format of tleap
        HBD.insert(pair<string,int>("LYS-HZ3",1)); //Additional line as the format of tleap
        HBD.insert(pair<string,int>("HIP-HD1",1)); //Additional line as the format of tleap
        HBD.insert(pair<string,int>("HIP-2HE",1));
        HBD.insert(pair<string,int>("HIP-HE2",1)); //Additional line as the format of tleap
//Neutral-N-HBD
        HBD.insert(pair<string,int>("LYN+2HZ",1));
        HBD.insert(pair<string,int>("LYN+3HZ",1));
        HBD.insert(pair<string,int>("LYN+HZ2",1)); //Additional line as the format of tleap
        HBD.insert(pair<string,int>("LYN+HZ3",1)); //Additional line as the format of tleap
        HBD.insert(pair<string,int>("ASN+1HD2",1));
        HBD.insert(pair<string,int>("ASN+2HD2",1));
        HBD.insert(pair<string,int>("GLN+1HE2",1));
        HBD.insert(pair<string,int>("GLN+2HE2",1));
        HBD.insert(pair<string,int>("ASN+HD21",1)); //Additional line as the format of tleap
        HBD.insert(pair<string,int>("ASN+HD22",1)); //Additional line as the format of tleap
        HBD.insert(pair<string,int>("GLN+HE21",1)); //Additional line as the format of tleap
        HBD.insert(pair<string,int>("GLN+HE22",1)); //Additional line as the format of tleap
//Neutral-O-HBD
        HBD.insert(pair<string,int>("THR*1HG",1));
        HBD.insert(pair<string,int>("THR*HG1",1));
        HBD.insert(pair<string,int>("SER*HG",1));
        HBD.insert(pair<string,int>("TYR*HH",1));
        HBD.insert(pair<string,int>("ASH*2HD",1));
        HBD.insert(pair<string,int>("ASH*HD2",1)); //Additional line as the format of tleap
        HBD.insert(pair<string,int>("GLH*2HE",1));
        HBD.insert(pair<string,int>("GLH*HE2",1)); //Additional line as the format of tleap
//Aromatic-HBD
        HBD.insert(pair<string,int>("TRP#HE1",1));
        HBD.insert(pair<string,int>("HIE#2HE",1));
        HBD.insert(pair<string,int>("HIE#HE2",1)); //Additional line as the format of tleap
        HBD.insert(pair<string,int>("HID#1HD",1));
        HBD.insert(pair<string,int>("HID#HD1",1)); //Additional line as the format of tleap

//Charged-HBA
        HBA.insert(pair<string,int>("ASP-OD1",1));
        HBA.insert(pair<string,int>("ASP-OD2",1));
        HBA.insert(pair<string,int>("GLU-OE1",1));
        HBA.insert(pair<string,int>("GLU-OE2",1));
//Neutral-N-HBA
        HBA.insert(pair<string,int>("LYN+NZ",1));
//Neutral-O-HBA
        HBA.insert(pair<string,int>("ASN*OD1",1));
        HBA.insert(pair<string,int>("GLN*OE1",1));
        HBA.insert(pair<string,int>("ASH*OD1",1));
        HBA.insert(pair<string,int>("ASH*OD2",1));
        HBA.insert(pair<string,int>("GLH*OE1",1));
        HBA.insert(pair<string,int>("GLH*OE2",1));
        HBA.insert(pair<string,int>("THR*OG1",1));
        HBA.insert(pair<string,int>("SER*OG",1));
//Aromatic-HBA
        HBA.insert(pair<string,int>("TYR#OH",1));
        HBA.insert(pair<string,int>("HIE#ND1",1));
        HBA.insert(pair<string,int>("HID#NE2",1));

        logP1.insert(pair<string,float>("ALA",1.0262));
        logP1.insert(pair<string,float>("VAL",1.6623));
        logP1.insert(pair<string,float>("LEU",2.0524));
        logP1.insert(pair<string,float>("ILE",2.0524));
        logP1.insert(pair<string,float>("PRO",0.3698));
        logP1.insert(pair<string,float>("MET",1.7594));
        logP2.insert(pair<string,float>("TRP",2.7303));
        logP2.insert(pair<string,float>("PHE",2.2490));
        logP2.insert(pair<string,float>("TYR",1.955));  
        logP2.insert(pair<string,float>("HID",0.972));
	logP2.insert(pair<string,float>("HIE",0.972));

        pmr1.insert(pair<string,float>("ALA",21.285));
        pmr1.insert(pair<string,float>("VAL",30.449));
        pmr1.insert(pair<string,float>("LEU",35.066));
        pmr1.insert(pair<string,float>("ILE",35.066));
        pmr1.insert(pair<string,float>("PRO",28.660));
        pmr1.insert(pair<string,float>("MET",38.610));
        pmr1.insert(pair<string,float>("GLY",16.690));
        pmr1.insert(pair<string,float>("CYS",29.464));
        pmr1.insert(pair<string,float>("THR",27.292));
        pmr1.insert(pair<string,float>("SER",22.697));
	pmr2.insert(pair<string,float>("HID",37.903));
	pmr2.insert(pair<string,float>("HIE",37.903));
        pmr2.insert(pair<string,float>("TRP",57.614));
        pmr2.insert(pair<string,float>("PHE",45.757));
        pmr2.insert(pair<string,float>("TYR",47.422));
}
