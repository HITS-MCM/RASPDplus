#!/bin/bash

#>  Copyright (c) 2020
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
#>  Manuscript in preparation, 2020

#>  Author of this script: Goutam Mukherjee
#>  Input information needed: (i) JobID, and (ii) user specified query scaffolds in a file named “scaffolds.txt”. 
#>  Please note that “scaffolds.txt” file must be present where the script, "scaffolds_search.sh" is executed.
#>  Purpose:
#>  Searching query scaffolds from a target SMILES file (target.smi). 

#>  RASPD+ screen million or customized small molecule databases against a target protein and final output is a file (FinalResult.txt) which contains predicted binding free energies of these small molecules and SMILES codes for small molecules (target.smi).

#>  These SMILES codes of samll molecule may contains several scaffolds/functional groups. If, one need to select a particular active scaffold from these databases, the SMILES code of this query active scaffolds  needs to  supply as a file name (scaffolds.txt). Please DON’T give file name other than scaffolds.txt.

#>  Command to run:
#>  bash scaffolds_search.sh <JobID>
#>  Say,
#>  bash scaffolds_search.sh 24190489_million_1NHZ_486_all
#>  This script may takes several minutes to hours depending on number of scaffolds present in the query file i.e.; in “scaffolds.txt” file.
#>  Output files are: (i) target_scaffold.smi; (ii) target_scaffold_be.txt
#>  Output file can be found inside the "JobID" folder

eval "$($conda_root/condabin/conda shell.bash hook)"

if [ -z "$raspd_root" ]
then
	echo -e '\033[1mPlease set path of the RASPD+ repository, miniconda and TRAPP to <path_to_RASPD+_repository>/config/init.sh file and source it\033[0m'
        echo -e '\033[1m   source <path_to_RASPD+_repository>/config/init.sh\033[0m'
        echo " "
        echo -e '\033[1m script usage:\033[0m'
        echo "bash scaffolds_search.sh <JobID>"
        exit
fi

export scripts_path=$raspd_root/scripts/chem
date

cwd=$1
if [ ! -d "$cwd" ]
  then
    echo -e '\033[7mPlease proide the full-path of the folder\033[0m'
    else if [ ! -f "scaffolds.txt" ]
        then
		echo -e '\033[1mPlease provide the input query scaffolds as a file name scaffolds.txt. Please keep the scaffolds.txt file at the same place where the scaffolds_search.sh is executed.\033[0m'
        else if [ ! -f "$cwd/target.smi" ]
	then
		echo -e '\033[1mPlease provide the target SMILES codes as a file name target.smi. Please keep the target.smi file inside the '$cwd/' folder\033[0m'
	else
	export download_path=`pwd`
##cwd=`find ~/ -name "$jobid"`
	cd $cwd	
        rm -rf scaffolds
        mkdir scaffolds
        cd scaffolds
	cp $download_path/scaffolds.txt .
	ln -s $scripts_path/substructure_search_rdkit.py . 
        seq 1 `wc -l scaffolds.txt | awk '{print($1)}'` >lst
        exec 3<scaffolds.txt
        exec 4<lst
                while read -r line <&3
                read -r lne <&4
        do
                echo $lne
                sed -n `echo $lne`p scaffolds.txt >$lne".smi"
	      conda activate raspdml
              python substructure_search_rdkit.py -q $lne".smi" -t ../target.smi -o $lne"_output.txt"
	      conda deactivate
              v=`wc -l $lne"_output.txt" | awk '{print($1)}'`
        zero=0;
        if [ $v -eq $zero ]; then
        echo "Error in input SMILES $lne" >$lne"_output.txt"
        fi
      if [ $v -ne $zero ]; then
              egrep -n "True" $lne"_output.txt" | awk -F \: '{print($1+1)}' >$lne".id"
              fi

      cat $lne".id" >>Scaffolds_Choice_by_User
        done
	sort -n Scaffolds_Choice_by_User | uniq >Scaffolds_Choice_by_User.txt
	cd ../
	
	if [ ! -f "FinalResult.txt" ]
then
	/usr/bin/perl $scripts_path/line_select.pl target.smi scaffolds/Scaffolds_Choice_by_User.txt >target_scaffold.smi
	echo -e "\e[1;31m Output file: $cwd/target_scaffold.smi"
	echo -e '\033[1m Please note here target_scaffold.smi, the output file, contains only thoes SMILES strings which have same substructure w.r.t. the query SMILES string, scaffolds.txt\033[0m'
else
	/usr/bin/perl $scripts_path/line_select.pl FinalResult.txt scaffolds/Scaffolds_Choice_by_User.txt >target_scaffold_be.txt
	echo -e "\e[1;31m Output file (i): $cwd/target_scaffold.smi"
	echo -e '\033[1m Please note here target_scaffold.smi, the first output file, contains only thoes SMILES strings which have same substructure w.r.t. the query SMILES string, scaffolds.txt\033[0m'
	echo -e "\e[1;31m Output file (ii): $cwd/target_scaffold_be.txt"
	echo -e '\033[1m Please note here target_scaffold_be.txt, the second output file, contains the predicted binding free energies of thoes SMILES strings which have same substructure w.r.t. the query SMILES string, scaffolds.txt\033[0m'
fi
date
fi
fi
fi
