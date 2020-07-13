#!/bin/bash

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
#>  How to run:
#>  bash single_molecule_scanning.sh <ligand.pdb> <Protein.pdb> <Identifier ID> <Method name>
#>  Say, lgand.pdb = lig.pdb
#>  Protein.pdb = 1NHZ.pdb
#>  Identifier ID =  486
#>  Method(s) = erf, rf, dnn, knn, svr, esvr, lr and a combinations of all of these methods (all).   
#>  bash single_molecule_scanning.sh lig.pdb 1NHZ.pdb 486 lr
#>  or,
#>  bash single_molecule_scanning.sh lig.pdb 1NHZ.pdb 486 erf

#>  Please note that all the input files MUST BE at the current working directory.
#>  Please note that input ligand size should be less than 300 atoms.
#>  All the jobs will be executed at the current working directory.
#>  Output file name is: FinalResult.txt
#>  Output file is there inside the JobID folder which is created while running this script.


if [ -z "$raspd_root" ]
then
	echo -e '\033[1mPlease set path of the RASPD+ repository, miniconda and TRAPP to <path_to_RASPD+_repository>/config/init.sh file and source it\033[0m'
        echo -e '\033[1m   source <path_to_RASPD+_repository>/config/init.sh\033[0m'
        echo " "
        echo -e '\033[1m script usage:\033[0m'
	echo "bash single_molecule_scanning.sh <ligand.pdb> <Protein.pdb> <Identifier ID> <Method name>"
        exit
fi

eval "$($conda_root/condabin/conda shell.bash hook)"
ligand=$1
protein=$2
hetid=$3
method=$4
min_prot_size=20;
if [  -f "$protein" ]
then
	prot_atom=`egrep "ATOM" $protein | egrep "CA" | egrep "ALA|ARG|ASH|ASN|ASP|CYM|CYS|CYX|GLH|GLN|GLU|GLH|GLY|HIS|HID|HIE|HIP|ILE|LEU|LYN|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL" | wc -l | awk '{print($1)}'`
fi

if [ -z "$ligand" ]
  then
    echo -e '\033[7mNo input ligand file ID is provided\033[0m'
        else if [ ! -f "$ligand" ]
        then
		echo -e '\033[1mInput ligand do not exist at the present working directory.\033[0m'
		echo -e '\033[1mPlease provide a valid input ligand-id with .pdb/.mol2/.sdf extension.\033[0m'
		echo -e '\033[1mAlternatively, if you have SMILES code instead of 3D coordinates of ligand molecule, then please provide it (SMILES strings) with a file name with a .txt extension\033[0m'
                        else if [ -z "$protein" ]
                                then
					echo -e '\033[7mNo input protein file ID is provided\033[0m'
                                                else if [ ! -f "$protein" ]
                                                        then
								echo -e '\033[1mPDB-ID of input protein do not exist at the present working directory\033[0m'
												else if [ $prot_atom -lt $min_prot_size ];
												then
													echo -e '\033[4mProtein contains less than 20 amino acid residues. Please provide a valid protein file\033[0m'
                                                                        				else if [ -z "$hetid" ]
                                                                                			then
														echo -e '\033[1mNo Identifier ID is provided\033[0m'
                                                                                                		else  if [ "$method" != 'erf' ] && [ "$method" != 'rf' ] && [ "$method" != 'dnn' ] && [ "$method" != 'knn' ] && [ "$method" != 'svr' ] && [ "$method" != 'esvr' ] && [ "$method" != 'lr' ] && [ "$method" != 'all' ]
                                                                                                        	then
															echo -e '\033[7mPlease choose a valid method name (e.g.; erf or, rf, or, dnn or, knn. or, svr, or, esvr, or, lr or, all)\033[0m'
                                                                                                                        	else
ml_scoring()
{
export ml_scripts_path=$raspd_root/scripts/ml
export ml_data_path=$raspd_root/weights

mkdir -p ml_runs
conda activate raspdml

if [ "$method" = "all" ]
then
        python $ml_scripts_path/infer.py -d ml_runs -w $ml_data_path -i Protein-Ligand_Parameter.txt -m erf rf dnn knn svr esvr lr -c 32
        conda deactivate
paste name ml_runs/erf.out ml_runs/rf.out ml_runs/dnn.out ml_runs/knn.out ml_runs/svr.out ml_runs/esvr.out ml_runs/lr.out >FinalResult
awk '{printf "%-20s;%0.3f;%0.3f;%0.3f;%0.3f;%0.3f;%0.3f;%0.3f\n", $1, $2, $4, $6, $8, $10, $12, $14}' FinalResult >FinalResult.txt
cat name_discard >temp
echo "molecule-id;erf;rf;dnn;knn;svr;esvr;lr" >>temp
cat FinalResult.txt >>temp
mv temp FinalResult.txt
fi

if [ "$method" = "erf" ]
then
        python $ml_scripts_path/infer.py -d ml_runs -w $ml_data_path -i Protein-Ligand_Parameter.txt -m erf -c 32
        conda deactivate
paste name ml_runs/erf.out >FinalResult
awk '{printf "%-20s;%0.3f\n", $1, $2}' FinalResult >FinalResult.txt
cat name_discard >temp
cat FinalResult.txt >>temp
mv temp FinalResult.txt
fi

if [ "$method" = "rf" ]
then
        python $ml_scripts_path/infer.py -d ml_runs -w $ml_data_path -i Protein-Ligand_Parameter.txt -m rf -c 32
        conda deactivate
paste name ml_runs/rf.out >FinalResult
awk '{printf "%-20s;%0.3f\n", $1, $2}' FinalResult >FinalResult.txt
cat name_discard >temp
cat FinalResult.txt >>temp
mv temp FinalResult.txt
fi

if [ "$method" = "dnn" ]
then
        python $ml_scripts_path/infer.py -d ml_runs -w $ml_data_path -i Protein-Ligand_Parameter.txt -m dnn -c 32
        conda deactivate
paste name ml_runs/dnn.out >FinalResult
awk '{printf "%-20s;%0.3f\n", $1, $2}' FinalResult >FinalResult.txt
cat name_discard >temp
cat FinalResult.txt >>temp
mv temp FinalResult.txt
fi

if [ "$method" = "knn" ]
then
        python $ml_scripts_path/infer.py -d ml_runs -w $ml_data_path -i Protein-Ligand_Parameter.txt -m knn -c 32
        conda deactivate
paste name ml_runs/knn.out >FinalResult
awk '{printf "%-20s;%0.3f\n", $1, $2}' FinalResult >FinalResult.txt
cat name_discard >temp
cat FinalResult.txt >>temp
mv temp FinalResult.txt
fi

if [ "$method" = "svr" ]
then
        python $ml_scripts_path/infer.py -d ml_runs -w $ml_data_path -i Protein-Ligand_Parameter.txt -m svr -c 32
        conda deactivate
paste name ml_runs/svr.out >FinalResult
awk '{printf "%-20s;%0.3f\n", $1, $2}' FinalResult >FinalResult.txt
cat name_discard >temp
cat FinalResult.txt >>temp
mv temp FinalResult.txt
fi

if [ "$method" = "esvr" ]
then
        python $ml_scripts_path/infer.py -d ml_runs -w $ml_data_path -i Protein-Ligand_Parameter.txt -m esvr -c 32
        conda deactivate
paste name ml_runs/esvr.out >FinalResult
awk '{printf "%-20s;%0.3f\n", $1, $2}' FinalResult >FinalResult.txt
cat name_discard >temp
cat FinalResult.txt >>temp
mv temp FinalResult.txt
fi

if [ "$method" = "lr" ]
then
        python $ml_scripts_path/infer.py -d ml_runs -w $ml_data_path -i Protein-Ligand_Parameter.txt -m lr -c 32
        conda deactivate
paste name ml_runs/lr.out >FinalResult
awk '{printf "%-20s;%0.3f\n", $1, $2}' FinalResult >FinalResult.txt
cat name_discard >temp
cat FinalResult.txt >>temp
mv temp FinalResult.txt
fi
}
cat FinalResult.txt >>../FinalResult.txt

export input_path=$raspd_root/bin
export download_path=`pwd`
export scripts_path=$raspd_root/scripts/chem
export data_path=$raspd_root/data

$input_path/pdb2chaininfo.exe $protein $hetid | sort -d | uniq >ligandinfo
	sz=`wc -l ligandinfo | awk '{print($1)}'`
	zero=0;
	if [ $sz -eq $zero ];
 	 then
	 echo -e '\033[1mPlease provide a valid Identifier-ID\033[0m'
    	 else
if 
        [[ $download_path/$ligand =~ \.mol2$ ]]; then
        ligand=`echo $ligand | cut -f1 -d'.'`
        cwd=`pwd`
	v=`awk -v min=1000 -v max=90000000 'BEGIN{srand(); print int(min+rand()*(max-min+1))}'`
        mkdir "_"$v
	lig_size=`wc -l "$ligand".mol2"" | awk '{print($1)}'`
        max_size=300;
        if [ $lig_size -gt $max_size ];
        then
        echo "Input ligand is too large to calculate the parameters"
	else
	conda activate trappenv	
        babel -imol2 $ligand".mol2" -opdb -O $ligand".pdb"
	conda deactivate
        $input_path/pdbtopdb.exe $ligand".pdb" "_"$v/$ligand".pdb"
        rm $ligand".pdb"
        cp $download_path/$protein "_"$v
        cd "_"$v
	echo $ligand >name
        $input_path/pdb2chaininfo.exe $protein $hetid | sort -d | uniq >ligandinfo
        egrep "`head -1 ligandinfo`" $protein | egrep -v "REMARK|END" >ligand.pdb
        egrep "ATOM" $protein | egrep "ALA|ARG|ASH|ASN|ASP|CYM|CYS|CYX|GLH|GLN|GLU|GLH|GLY|HIS|HID|HIE|HIP|ILE|LEU|LYN|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL" >complex
	$input_path/pdb2noh.exe complex >complex_noh.pdb
	conda activate trappenv
	tleap -f $data_path/s1.cmd >leap.log
	conda deactivate
	egrep -v "REMARK|END" complex >complex.pdb
        cp complex.pdb protein.pdb
        $input_path/pdb2pdb.exe ligand.pdb >>complex.pdb  
        $input_path/CM-PL.exe complex.pdb  >com_complex
        $input_path/CM-PL.exe $ligand".pdb" >com_ligand
        $input_path/maxD_calculator.exe $ligand".pdb" com_ligand | sort -n | tail -n 1 >maxdistance
        paste maxdistance maxdistance >DMax
        $input_path/protein-res-cm-cal.exe protein.pdb >protein.ent
        $input_path/protein_param_calc.exe protein.ent DMax `cat com_complex` >>Protein_Parameter
        $scripts_path/trapp_volume.sh complex >>leap.log
	mkdir LigParm
	cp $ligand".pdb" LigParm/
	cd LigParm/
        $scripts_path/RASPD_Ligand_parameter.sh $ligand
	cd ../
        paste Protein_Parameter volm LigParm/MOLECULARPROPERTY LigParm/Molecular_Weight >Protein-Ligand_Parameter
	sed -i 's/inf/0/g' Protein-Ligand_Parameter
        $input_path/paste_ml_format.exe Protein-Ligand_Parameter >Protein-Ligand_Parameter.txt
        var=`egrep "0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0" Protein-Ligand_Parameter.txt | wc -l`
        zero=0;
        if [ $var -eq $zero ]; then
        ml_scoring       ##Output is FinalResult.txt
        else
	echo "Structural parameters of input $ligand".mol2" is not optimized. Please optimized the ligand molecule and resubmit the job." >FinalResult.txt
        fi
	mkdir tmpdir
	mv * tmpdir/
	cp tmpdir/FinalResult.txt .
	cd ../
	prot=`echo $protein | cut -f1 -d'.'`
        mv "_"$v $v"_mol2_screening_"$prot"_"$hetid"_"$method
	echo -e "\e[1;31m Output directory: $v"_mol2_screening_"$prot"_"$hetid"_"$method"
fi fi

if 
        [[ $download_path/$ligand =~ \.sdf$ ]]; then
        ligand=`echo $ligand | cut -f1 -d'.'`
        cwd=`pwd`
        v=`awk -v min=1000 -v max=90000000 'BEGIN{srand(); print int(min+rand()*(max-min+1))}'`
        mkdir "_"$v       
        lig_size=`wc -l "$ligand.sdf" | awk '{print($1)}'`
        max_size=300;
        if [ $lig_size -gt $max_size ];
        then
        echo "Input ligand is too large to calculate the parameters"
        else	
	conda activate trappenv
	babel -isdf $ligand".sdf" -opdb -O $ligand".pdb"
        $input_path/pdbtopdb.exe $ligand".pdb" "_"$v/$ligand".pdb"
        rm $ligand".pdb"
        cp $download_path/$protein "_"$v/
        cd "_"$v/
	echo $ligand >name
        $input_path/pdb2chaininfo.exe $protein $hetid | sort -d | uniq >ligandinfo
        egrep "`head -1 ligandinfo`" $protein | egrep -v "REMARK|END" >ligand.pdb
        egrep "ATOM" $protein | egrep "ALA|ARG|ASH|ASN|ASP|CYM|CYS|CYX|GLH|GLN|GLU|GLH|GLY|HIS|HID|HIE|HIP|ILE|LEU|LYN|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL" >complex
	$input_path/pdb2noh.exe complex >complex_noh.pdb
        conda activate trappenv
        tleap -f $data_path/s1.cmd >leap.log
        conda deactivate	
	egrep -v "REMARK|END" complex >complex.pdb
        cp complex.pdb protein.pdb
        $input_path/pdb2pdb.exe ligand.pdb >>complex.pdb  
        $input_path/CM-PL.exe complex.pdb  >com_complex
        $input_path/CM-PL.exe $ligand".pdb" >com_ligand
        $input_path/maxD_calculator.exe $ligand".pdb" com_ligand | sort -n | tail -n 1 >maxdistance
        paste maxdistance maxdistance >DMax
        $input_path/protein-res-cm-cal.exe protein.pdb >protein.ent
        $input_path/protein_param_calc.exe protein.ent DMax `cat com_complex` >>Protein_Parameter
        $scripts_path/trapp_volume.sh complex >>leap.log
        mkdir LigParm
        cp $ligand".pdb" LigParm/
        cd LigParm/
        $scripts_path/RASPD_Ligand_parameter.sh $ligand
        cd ../
        paste Protein_Parameter volm LigParm/MOLECULARPROPERTY LigParm/Molecular_Weight >Protein-Ligand_Parameter
	sed -i 's/inf/0/g' Protein-Ligand_Parameter
	$input_path/paste_ml_format.exe Protein-Ligand_Parameter >Protein-Ligand_Parameter.txt
        var=`egrep "0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0" Protein-Ligand_Parameter.txt | wc -l`
        zero=0;
        if [ $var -eq $zero ]; then
        ml_scoring               ##Output is FinalResult.txt
        else
        echo "Structural parameters of input $ligand".sdf" is not optimized. Please optimized the ligand molecule and resubmit the job." >FinalResult.txt
        fi
	mkdir tmpdir
        mv * tmpdir/
        cp tmpdir/FinalResult.txt .
	cd ../
	prot=`echo $protein | cut -f1 -d'.'`
        mv "_"$v $v"_sdf_screening_"$prot"_"$hetid"_"$method
	echo -e "\e[1;31m Output directory: $v"_sdf_screening_"$prot"_"$hetid"_"$method"
fi fi

if 
        [[ $download_path/$ligand =~ \.pdb$ ]]; then
        ligand=`echo $ligand | cut -f1 -d'.'`
        cwd=`pwd`
        v=`awk -v min=1000 -v max=90000000 'BEGIN{srand(); print int(min+rand()*(max-min+1))}'`
        mkdir "_"$v
        $input_path/pdbtopdb.exe $download_path/$ligand".pdb" "_"$v/$ligand".pdb"
        cp $download_path/$protein "_"$v/
        cd "_"$v
	lig_size=`wc -l "$ligand.pdb" | awk '{print($1)}'`
        max_size=300;
        if [ $lig_size -gt $max_size ];
        then
	echo "Input ligand is too large to calculate the parameters"
	else
	echo $ligand >name
        $input_path/pdb2chaininfo.exe $protein $hetid | sort -d | uniq >ligandinfo
        egrep "`head -1 ligandinfo`" $protein | egrep -v "REMARK|END" >ligand.pdb
        egrep "ATOM" $protein | egrep "ALA|ARG|ASH|ASN|ASP|CYM|CYS|CYX|GLH|GLN|GLU|GLH|GLY|HIS|HID|HIE|HIP|ILE|LEU|LYN|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL" >complex
	$input_path/pdb2noh.exe complex >complex_noh.pdb
        conda activate trappenv
        tleap -f $data_path/s1.cmd >leap.log
        conda deactivate        
	egrep -v "REMARK|END" complex >complex.pdb
        cp complex.pdb protein.pdb
        $input_path/pdb2pdb.exe ligand.pdb >>complex.pdb  
        $input_path/CM-PL.exe complex.pdb  >com_complex
        $input_path/CM-PL.exe $ligand".pdb" >com_ligand
        $input_path/maxD_calculator.exe $ligand".pdb" com_ligand | sort -n | tail -n 1 >maxdistance
        paste maxdistance maxdistance >DMax
        $input_path/protein-res-cm-cal.exe protein.pdb >protein.ent
        $input_path/protein_param_calc.exe protein.ent DMax `cat com_complex` >>Protein_Parameter
        $scripts_path/trapp_volume.sh complex >>leap.log
        mkdir LigParm
        cp $ligand".pdb" LigParm/
        cd LigParm/
        $scripts_path/RASPD_Ligand_parameter.sh $ligand
        cd ../
        paste Protein_Parameter volm LigParm/MOLECULARPROPERTY LigParm/Molecular_Weight >Protein-Ligand_Parameter
	sed -i 's/inf/0/g' Protein-Ligand_Parameter
	$input_path/paste_ml_format.exe Protein-Ligand_Parameter >Protein-Ligand_Parameter.txt
        var=`egrep "0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0" Protein-Ligand_Parameter.txt | wc -l`
        zero=0;
        if [ $var -eq $zero ]; then
        ml_scoring          ##Output is FinalResult.txt
        else
        echo "Structural parameters of input $ligand".pdb" is not optimized. Please optimized the ligand molecule and resubmit the job." >FinalResult.txt
        fi
	mkdir tmpdir
        mv * tmpdir/
        cp tmpdir/FinalResult.txt .
	cd ../
	prot=`echo $protein | cut -f1 -d'.'`
	mv "_"$v $v"_pdb_screening_"$prot"_"$hetid"_"$method
	echo -e "\e[1;31m Output directory: $v"_pdb_screening_"$prot"_"$hetid"_"$method"
fi fi

if 
        [[ $download_path/$ligand =~ \.txt$ ]]; then
        ligand=`echo $ligand | cut -f1 -d'.'`
        echo $ligand
	lig_size=`wc -l "$ligand.txt" | awk '{print($1)}'`
        size=0;
	if [ $lig_size -eq $size ]; then
		echo "Please provide atleast one SMILES code in the input scaffolds.txt file"
		else
        cwd=`pwd`
        v=`awk -v min=1000 -v max=90000000 'BEGIN{srand(); print int(min+rand()*(max-min+1))}'`
        mkdir "_"$v
	cp $download_path/$ligand".txt" "_"$v/$ligand".txt"
	cp $download_path/$protein "_"$v/
	cd "_"$v
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
	echo -e '\033[1mPlease ignore the Segmentation fault      (core dumped) error if any\033[0m'
	echo -e '\033[1mThis error may appeared if the output of Connect2.0.exe can not read by the windex_single.exe, a programme to calculate wiener index of a small molecule or other programme\033[0m'
	echo -e '\033[1mYou may get this error due to the problem in the input structural parameters as well\033[0m'
	echo " "

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
	paste LigParam/maxdistance LigParam/maxdistance >DMax
        $input_path/pdb2chaininfo.exe $protein $hetid | sort -d | uniq >ligandinfo
        egrep "`head -1 ligandinfo`" $protein | egrep -v "REMARK|END" >ligand.pdb
        egrep "ATOM" $protein | egrep "ALA|ARG|ASH|ASN|ASP|CYM|CYS|CYX|GLH|GLN|GLU|GLH|GLY|HIS|HID|HIE|HIP|ILE|LEU|LYN|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL" >complex
	$input_path/pdb2noh.exe complex >complex_noh.pdb
        conda activate trappenv
        tleap -f $data_path/s1.cmd >leap.log
        conda deactivate
	egrep -v "REMARK|END" complex >complex.pdb
        cp complex.pdb protein.pdb
        $input_path/pdb2pdb.exe ligand.pdb >>complex.pdb  
        $input_path/CM-PL.exe complex.pdb  >com_complex	
        $input_path/protein-res-cm-cal.exe protein.pdb >protein.ent
        $input_path/protein_param_calc.exe protein.ent DMax `cat com_complex` >>Protein_Parameter
        $scripts_path/trapp_volume.sh complex >>leap.log
	$input_path/print_volm_multiple.exe `wc -l smiles.corr | awk '{print($1)}'` `cat volm` >VOLUME
        paste smiles.corr Protein_Parameter VOLUME LigParam/MOLECULARPROPERTY LigParam/Molecular_Weight >Protein-Ligand_Parameter
	sed -i 's/inf/0/g' Protein-Ligand_Parameter
        $input_path/paste_ml_format_4smiles.exe Protein-Ligand_Parameter >Protein-Ligand_Parameter.txt
        var=`egrep "0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0" Protein-Ligand_Parameter.txt | wc -l`
        count=`wc -l lst | awk '{print($1)}'`;
        zero=0;
	if [ $var -eq $count ]; then
        echo "Error in input SMILES codes" >FinalResult.txt
        fi

	if [ $var -eq $zero ]; then
        sed -n 2,1000000000p Protein-Ligand_Parameter.txt | awk -F \; '{print($1)}' >name
	ml_scoring
	fi

	if [ $var -lt $count ]; then
	egrep -v "0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0" Protein-Ligand_Parameter.txt >Protein-Ligand_Parameter_nonzero.txt
	if [ -f smiles.err ]; then
	awk '{print("SMILES-ID:" $1 "is discarded. Please check the input SMILES string.")}' smiles.err >name_discard
	fi
	egrep -n "0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0" Protein-Ligand_Parameter.txt | awk -F":" '{print("SMILES-ID:" $1-1" is discarded due to the problem in the structural parameters.")}' >>name_discard
	mv Protein-Ligand_Parameter_nonzero.txt Protein-Ligand_Parameter.txt
	sed -n 2,1000000000p Protein-Ligand_Parameter.txt | awk -F \; '{print($1)}' >name
	ml_scoring
	fi
	mkdir tmpdir
        mv * tmpdir/
        cp tmpdir/FinalResult.txt .
	cd ../
	prot=`echo $protein | cut -f1 -d'.'`
        mv "_"$v $v"_smile_screening_"$prot"_"$hetid"_"$method
	echo -e "\e[1;31m Output directory: $v"_smile_screening_"$prot"_"$hetid"_"$method"
fi fi
									fi 
								fi 
							fi 
						fi 	
					fi 
				fi 
			fi 
		fi
