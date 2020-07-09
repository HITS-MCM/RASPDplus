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
#>  Manuscript in preparation, 2020

#>  Author of this script: Goutam Mukherjee
#>  Purpose: Generate physicochemical parameters for customized small molecules. Parameters will be saved automatically in "$raspd_root/customized_data/" directory
#>  Input file is *.txt file. Say, molecule.txt.
#>  molecule.txt may contains one or more SMILES strings of ligand (drug-like) molecules
#>  How to run:
#>  bash lig_parameters_gen.sh molecule.txt <Yes (please note, file extension (*.txt) is mandatory)
#>  This *.txt file MUST BE present at the current working directory where, the job will be executed by a script name "lig_parameters_gen.sh".

if [ -z "$raspd_root" ]
then
	echo -e '\033[1mPlease set path of the RASPD+ repository, miniconda and TRAPP to <path_to_RASPD+_repository>/config/init.sh file and source it\033[0m'
        echo -e '\033[1m   source <path_to_RASPD+_repository>/config/init.sh\033[0m'
        echo " "
        echo -e '\033[1m script usage:\033[0m'
	echo "bash lig_parameters_gen.sh <filename.txt>"
        exit
fi
	echo -e '\033[1m This script will autoreplace the all parameter files which are saved into the <path_to_RASPDplus_repository>/customized_data/ folder\033[0m'
	echo -e '\033[1m It is strongly recommend before you run this script, please take a backup of all the files which were saved in <path_to_RASPDplus_repository>/customized_data/ folder\033[0m'
	echo -e '\033[7m Proceed to run?\033[0m'
read answer
	if [ "$answer" != "${answer#[Yy]}" ]
		then
			eval "$($conda_root/condabin/conda shell.bash hook)"

			ligand=$1
				if [ -z "$ligand" ]
				then
				echo -e '\033[7mNo input file ID is provided\033[0m'
        				else if [ ! -f "$ligand" ]
   				     	then
        			        echo -e '\033[1mInput SMILES of ligand molecules do not exist at the present working directory. Please provide it (SMILES strings) with a .txt extension\033[0m'
        			        	else
						export input_path=$raspd_root/bin
						export download_path=`pwd`
						export scripts_path=$raspd_root/scripts/chem
						export customized_data_path=$raspd_root/customized_data
							if 
						        [[ $download_path/$ligand =~ \.txt$ ]]; then
						        ligand=`echo $ligand | cut -f1 -d'.'`
						        echo $ligand
						        cwd=`pwd`
						        v=`awk -v min=1000 -v max=90000000 'BEGIN{srand(); print int(min+rand()*(max-min+1))}'`
						        mkdir "@"$v
							cp $download_path/$ligand".txt" "@"$v/$ligand".txt"
							cd "@"$v
							cp $ligand".txt" list
							seq 1 `wc -l $ligand".txt" | awk '{print($1)}'` >lst
							exec 3<list
							exec 4<lst
							while read -r line <&3
						        read -r lne <&4
							        do
								echo $lne
								sed -n `echo $lne`p list >$lne".smi"
								done
								echo -e '\033[1mConversion of SMILES to 3D structures by Open Babel (OpenBabel)\033[0m'
							exec 3<lst
							while read -r lne <&3
								do
								conda activate trappenv
								obabel -ismi $lne".smi" -opdb -O $lne".pdb" --gen3d
								var=`wc -l $lne".pdb" | awk '{print($1)}'`
								zero=0;
								if [ $var -eq $zero ]; then
									echo "$lne: No SMILES string" >>smiles.err
								else
									echo $lne >>smiles.corr
								fi
							done
							mkdir LigParam
							cd LigParam
							echo -e '\033[1mCalculations of physicochemical parameters\033[0m'
							exec 3<../smiles.corr
					                while read -r lne <&3
								do
								cp ../$lne".pdb" .
								egrep "ATOM|HETATM"  $lne".pdb" >tmp
								$input_path/pdbtopdb.exe tmp $lne".pdb"
						        	$input_path/CM-PL.exe $lne".pdb" >com_ligand
						        	$input_path/maxD_calculator.exe $lne".pdb" com_ligand | sort -n | tail -n 1 >>maxdistance
								$scripts_path/RASPD_Ligand_parameter.sh $lne			
							done
							cd ../

							paste smiles.corr LigParam/MOLECULARPROPERTY LigParam/Molecular_Weight LigParam/maxdistance >list_PARAMETER_MW_maxd_Customized	
						        sed -i 's/inf/0/g' list_PARAMETER_MW_maxd_Customized
      							$input_path/param_chek_nonzero.exe list_PARAMETER_MW_maxd_Customized $customized_data_path/PARAMETER_MW_Customized $customized_data_path/maxdistance_customized >$customized_data_path/list_customized
							cp list $customized_data_path/target_customized.smi
						        echo -e "\e[1;31m Calculations are over and parameters are autosaved at <path_to_RASPD+_repository>/customized_data/ folder"
							fi	
					fi 		
				fi
			else
			echo "Thanks!"
			fi

