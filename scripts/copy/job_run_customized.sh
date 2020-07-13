#>  Copyright (c) 2013 2014 2015 2016 2017 2018 2019 2020
#>  Heidelberg University and Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org)
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

#>  Author of this script: Goutam Mukherjee
#>  How to run this script:
#>  bash job_run_customized.sh 1NHZ 486 erf
#>  Here, 
#>  1NHZ is the PDBid
#>  486 is the Active_Site_Identifier of PDBid 1NHZ
#>  erf is the Method used

#>  Methods available for RASPD+ Screening: 
#>  Extremely Random Forest (erf) 
#>  Random Forest (rf) 
#>  Deep Neural Network (dnn) 
#>  k-Nearest neighbors (knn)
#>  linear Support Vector Regression (svr) 
#>  (non-linear) Epsilon Support Vector Regression (esvr) 
#>  Linear Regression (lr)
#>  A combinations of all the seven methods (all)
#>  Using this script it is possible to screen customized small molecules dataset against a protein.
#>  Output files: 
#>    (i) FinalResult.txt (Contains molecule-id in first column and predicted binding free energies of the customized data set in second column)
#>    (ii) target.smi (Contains SMILE Code of the customized data set)
#>  Output file is there inside the JobID folder to be created while running this script.

if [ -z "$raspd_root" ]
then
	echo -e '\033[1mPlease set path of the RASPD+ repository, miniconda and TRAPP to <path_to_RASPD+_repository>/config/init.sh file and source it\033[0m'
        echo -e '\033[1m   source <path_to_RASPD+_repository>/config/init.sh\033[0m'
	echo " "
        echo -e '\033[1m script usage:\033[0m'
        echo "bash job_run_customized.sh <protein file name> <Active_Site_Identifier> <Method>"
        exit
fi

export scripts_path=$raspd_root/scripts/chem
export data_path=$raspd_root/data

PDBid=$1
Active_Site_Identifier=$2
Method=$3
min_prot_size=20;
if [  -f "$PDBid".pdb"" ]
then
        prot_atom=`egrep "ATOM" $PDBid".pdb" | egrep "ALA|ARG|ASH|ASN|ASP|CYM|CYS|CYX|GLH|GLN|GLU|GLH|GLY|HIS|HID|HIE|HIP|ILE|LEU|LYN|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL" | egrep "CB" | wc -l | awk '{print($1)}'`
fi

v=`awk -v min=1000 -v max=90000000 'BEGIN{srand(); print int(min+rand()*(max-min+1))}'`

if [ -z "$PDBid" ]
  then
   echo -e '\033[7mNo input protein file ID is provided\033[0m'
        else if [ ! -f "$PDBid".pdb"" ]
        then
		echo -e '\033[1mPDB-ID of input protein do not exist at the present working directory\033[0m'
			else if [ $prot_atom -lt $min_prot_size ];
                             then
			     echo -e '\033[4mProtein contains less than 20 amino acid residues. Please provide a valid protein file\033[0m'
                		else if [ -z "$Active_Site_Identifier" ]
                        		then
					echo -e '\033[1mNo Identifier ID is provided\033[0m'
                                 		else  if [ "$Method" != 'erf' ] && [ "$Method" != 'rf' ] && [ "$Method" != 'dnn' ] && [ "$Method" != 'knn' ] && [ "$Method" != 'svr' ] && [ "$Method" != 'esvr' ] && [ "$Method" != 'lr' ] && [ "$Method" != 'all' ]
                                        	then
							echo -e '\033[7mPlease choose a valid method name (e.g.; erf or, rf, or, dnn or, knn. or, svr, or, esvr, or, lr or, all)\033[0m'
                                                	else
                                                        $scripts_path/script_customized.sh $v $PDBid $Active_Site_Identifier $Method `cat $data_path/select_parameter.txt | awk '{printf($2" ")}'`  ##[For screening customized small molecules dataset]
								egrep -v "HEADER|TITLE|COMPND|SOURCE|KEYWDS|AUTHOR|REMARK|JRNL|SEQADV|CONECT" $PDBid".pdb" | egrep "ATOM|HETATM" | egrep -v "ALA|ARG|ASH|ASN|ASP|CYM|CYS|CYX|GLH|GLN|GLU|GLH|GLY|HIS|HID|HIE|HIP|ILE|LEU|LYN|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|HOH|WAT|SOL" | awk '{print $4}' | sort -d | uniq >resid
								sz=`egrep "$Active_Site_Identifier" resid | wc -l`
								if [ $sz -ne 0 ]
							  	then
                                                        	mv "Customized_Screening_"$v $v"_customized_"$PDBid"_"$Active_Site_Identifier"_"$Method
								echo -e "\e[1;31m Output directory: "$v"_customized_"$PDBid"_"$Active_Site_Identifier"_"$Method
								rm resid
							else
								rm -rf "Customized_Screening_"$v
							fi
						fi
                                        fi
                                fi
                        fi
                fi
