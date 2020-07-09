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
#>  Manuscript in preparation, 2020

#>  Author of this script: Stefan Holderbach, Heidelberg University and Goutam Mukherjee
#>  How to run:
#>  bash script_customized.sh <JobID> <protein-4-letter-code (without “.pdb” extension)> <Active_Site_Identifier ID> <Method> <`cat $data_path/select_parameter.txt | awk '{printf($2" ")}'`>
#>  Say, JobID is: Customized_Screening_35955430 
#>  protein-4-letter-code = 1NHZ.pdb 
#>  Identifier ID =  486
#>  Method(s) = erf, rf, dnn, knn, svr, esvr, lr and a combinations of all of these methods (all).   
#>  bash script_customized.sh Customized_Screening_35955430 1NHZ 486 lr `cat $data_path/select_parameter.txt | awk '{printf($2" ")}'`
#>  or,
#>  bash script_customized.sh Customized_Screening_35955430 1NHZ 486 erf `cat $data_path/select_parameter.txt | awk '{printf($2" ")}'`
#>  Output files: 
#>    (i) FinalResult.txt (Contains predicted binding free energies of the million molecules)
#    (ii) target.smi (Contains SMILES Code of the million molecules)

date
eval "$($conda_root/condabin/conda shell.bash hook)"
jobid=`eval echo "\${1}"`
pdb=`eval echo "\${2}"`
hetid=`eval echo "\${3}"`
method=`eval echo "\${4}"`   		##Newly added
cores=`eval echo "\${5}"`
wimin=`eval echo "\${6}"`
wimax=`eval echo "\${7}"`
hbdmin=`eval echo "\${8}"`
hbdmax=`eval echo "\${9}"`
hbamin=`eval echo "\${10}"`
hbamax=`eval echo "\${11}"`
logpmin=`eval echo "\${12}"`
logpmax=`eval echo "\${13}"`
mrmin=`eval echo "\${14}"`
mrmax=`eval echo "\${15}"`
mwmin=`eval echo "\${16}"`
mwmax=`eval echo "\${17}"`
aff=`eval echo "\${18}"`

if [ -z "$raspd_root" ]
then
	echo -e '\033[1mPlease set path of the RASPD+ repository, miniconda and TRAPP to <path_to_RASPD+_repository>/config/init.sh file and source it\033[0m'
        echo -e '\033[1m   source <path_to_RASPD+_repository>/config/init.sh\033[0m'
        echo " "
        echo -e '\033[1m You should not run this script directly\033[0m'
        echo -e '\033[1m Please go to the ../copy/ folder and run the script job_run_customized.sh\033[0m'
        echo -e '\033[1m job_run_customized.sh call this script to calculate binding energy of customize molecules against a target protein\033[0m'
        exit
fi

export input_path=$raspd_root/bin
export download_path=`pwd`
export scripts_path=$raspd_root/scripts/chem
export customized_data_path=$raspd_root/customized_data
export ml_scripts_path=$raspd_root/scripts/ml
export ml_data_path=$raspd_root/weights


cwd=`pwd`		   											## Added on 14.06.19
mkdir "Customized_Screening_"$jobid		   											## Added on 14.06.19
cd $cwd/"Customized_Screening_"$jobid		   											## Added on 14.06.19

ln -s $input_path/* .
ln -s $scripts_path/* .
ln -s $download_path/$pdb".pdb" .                                                                               ## Added on 14.06.19
ln -s $customized_data_path . 

./pdb2chaininfo.exe $pdb".pdb" $hetid | sort -d | uniq >ligandinfo			##Select the Small Molecules ID based on user input information (hetid)
sz=`wc -l ligandinfo | awk '{print($1)}'`
zero=0;
if [ $sz -eq $zero ];
  then
	  echo -e '\033[1mPlease provide a valid Identifier-ID.\nAs per RASPD+ nomenclature, a Identifier-ID is a three letter code in uppercase other than the amino acid residues such as "ALA", "GLY" etc. or, water such as, "HOH", "WAT", "SOL".\nSay, Identifier-ID of 1NHZ (https://www.rcsb.org/structure/1NHZ) is 486\033[0m'
    else
egrep "`head -1 ligandinfo`" $pdb".pdb" | egrep -v "REMARK|END" >ligand.pdb			## Added on 14.06.19; select the coordinates of ligand molecule from the uploaded protein-ligand complex
egrep "ATOM" $pdb".pdb" | egrep "ALA|ARG|ASH|ASN|ASP|CYM|CYS|CYX|GLH|GLN|GLU|GLH|GLY|HIS|HID|HIE|HIP|ILE|LEU|LYN|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL" >complex
./pdb2noh.exe complex >complex_noh.pdb

conda activate trappenv
tleap -f $customized_data_path/s1.cmd >leap.log										## Added on 14.06.19; Need to change
conda deactivate
egrep -v "REMARK|END" complex >complex.pdb
./pdb2pdb.exe ligand.pdb >>complex.pdb							##Remove chain information from ligand and paste below the protein coordinates.
mv $pdb".pdb" $pdb".ent"
cp complex.pdb $pdb".pdb"

./trapp_volume.sh $pdb >>leap.log
./print_volm_customized.exe `wc -l $customized_data_path/maxdistance_customized | awk '{print($1)}'` `cat volm` >Volume_1M

##EndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEnd

## Million Molecules Scoring Calculation
./CM-PL.exe $pdb".pdb" >com
./protein-res-cm-cal.exe $pdb".pdb" | egrep -v "DRG" >protein.ent
./protein_param_calc.exe protein.ent $customized_data_path/maxdistance_customized `cat com` >>Protein_Parameter
paste Protein_Parameter Volume_1M $customized_data_path/PARAMETER_MW_Customized >PARAMETERS_MILLION
./paste.exe PARAMETERS_MILLION >PARAMETERS_MILLION.txt
conda activate raspdml
mkdir -p ml_runs

if [ "$method" = "all" ]
then
        python $ml_scripts_path/infer.py -d ml_runs -w $ml_data_path -i PARAMETERS_MILLION.txt -m erf rf dnn knn svr esvr lr -c $cores 
        conda deactivate
paste ml_runs/erf.out ml_runs/rf.out ml_runs/dnn.out ml_runs/knn.out ml_runs/svr.out ml_runs/esvr.out ml_runs/lr.out $customized_data_path/PARAMETER_MW_Customized $customized_data_path/list_customized >pred_be_all
./final_select_all.exe pred_be_all $wimin $wimax $hbdmin $hbdmax $hbamin $hbamax $logpmin $logpmax $mrmin $mrmax $mwmin $mwmax $aff PDBID >FinalResult.txt
perl line_select.pl $customized_data_path/target_customized.smi PDBID >target.smi
fi

if [ "$method" = "erf" ]
then
        python $ml_scripts_path/infer.py -d ml_runs -w $ml_data_path -i PARAMETERS_MILLION.txt -m erf -c $cores 
        conda deactivate
paste ml_runs/erf.out $customized_data_path/PARAMETER_MW_Customized $customized_data_path/list_customized >pred_be_erf
./final_select.exe pred_be_erf $wimin $wimax $hbdmin $hbdmax $hbamin $hbamax $logpmin $logpmax $mrmin $mrmax $mwmin $mwmax $aff PDBID >FinalResult.txt
##sort -t ";" -k2n Result | awk '{print($1$2$3)}' >FinalResult.txt
perl line_select.pl $customized_data_path/target_customized.smi PDBID >target.smi
fi

if [ "$method" = "rf" ]
then
        python $ml_scripts_path/infer.py -d ml_runs -w $ml_data_path -i PARAMETERS_MILLION.txt -m rf -c $cores
        conda deactivate
paste ml_runs/rf.out $customized_data_path/PARAMETER_MW_Customized $customized_data_path/list_customized >pred_be_rf
./final_select.exe pred_be_rf $wimin $wimax $hbdmin $hbdmax $hbamin $hbamax $logpmin $logpmax $mrmin $mrmax $mwmin $mwmax $aff PDBID >FinalResult.txt
##sort -t ";" -k2n Result | awk '{print($1$2$3)}' >FinalResult.txt
perl line_select.pl $customized_data_path/target_customized.smi PDBID >target.smi
fi

if [ "$method" = "dnn" ]
then
        python $ml_scripts_path/infer.py -d ml_runs -w $ml_data_path -i PARAMETERS_MILLION.txt -m dnn -c $cores
        conda deactivate
paste ml_runs/dnn.out $customized_data_path/PARAMETER_MW_Customized $customized_data_path/list_customized >pred_be_dnn
./final_select.exe pred_be_dnn $wimin $wimax $hbdmin $hbdmax $hbamin $hbamax $logpmin $logpmax $mrmin $mrmax $mwmin $mwmax $aff PDBID >FinalResult.txt
##sort -t ";" -k2n Result | awk '{print($1$2$3)}' >FinalResult.txt
perl line_select.pl $customized_data_path/target_customized.smi PDBID >target.smi
fi

if [ "$method" = "knn" ]
then
        python $ml_scripts_path/infer.py -d ml_runs -w $ml_data_path -i PARAMETERS_MILLION.txt -m knn -c $cores
        conda deactivate
paste ml_runs/knn.out $customized_data_path/PARAMETER_MW_Customized $customized_data_path/list_customized >pred_be_knn
./final_select.exe pred_be_knn $wimin $wimax $hbdmin $hbdmax $hbamin $hbamax $logpmin $logpmax $mrmin $mrmax $mwmin $mwmax $aff PDBID >FinalResult.txt
##sort -t ";" -k2n Result | awk '{print($1$2$3)}' >FinalResult.txt
perl line_select.pl $customized_data_path/target_customized.smi PDBID >target.smi
fi

if [ "$method" = "svr" ]
then
        python $ml_scripts_path/infer.py -d ml_runs -w $ml_data_path -i PARAMETERS_MILLION.txt -m svr -c $cores
        conda deactivate
paste ml_runs/svr.out $customized_data_path/PARAMETER_MW_Customized $customized_data_path/list_customized >pred_be_svr
./final_select.exe pred_be_svr $wimin $wimax $hbdmin $hbdmax $hbamin $hbamax $logpmin $logpmax $mrmin $mrmax $mwmin $mwmax $aff PDBID >FinalResult.txt
##sort -t ";" -k2n Result | awk '{print($1$2$3)}' >FinalResult.txt
perl line_select.pl $customized_data_path/target_customized.smi PDBID >target.smi
fi

if [ "$method" = "esvr" ]
then
        python $ml_scripts_path/infer.py -d ml_runs -w $ml_data_path -i PARAMETERS_MILLION.txt -m esvr -c $cores
        conda deactivate
paste ml_runs/esvr.out $customized_data_path/PARAMETER_MW_Customized $customized_data_path/list_customized >pred_be_esvr
./final_select.exe pred_be_esvr $wimin $wimax $hbdmin $hbdmax $hbamin $hbamax $logpmin $logpmax $mrmin $mrmax $mwmin $mwmax $aff PDBID >FinalResult.txt
##sort -t ";" -k2n Result | awk '{print($1$2$3)}' >FinalResult.txt
perl line_select.pl $customized_data_path/target_customized.smi PDBID >target.smi
fi

if [ "$method" = "lr" ]
then
        python $ml_scripts_path/infer.py -d ml_runs -w $ml_data_path -i PARAMETERS_MILLION.txt -m lr -c $cores
        conda deactivate
paste ml_runs/lr.out $customized_data_path/PARAMETER_MW_Customized $customized_data_path/list_customized >pred_be_lr
./final_select.exe pred_be_lr $wimin $wimax $hbdmin $hbdmax $hbamin $hbamax $logpmin $logpmax $mrmin $mrmax $mwmin $mwmax $aff PDBID >FinalResult.txt
##sort -t ";" -k2n Result | awk '{print($1$2$3)}' >FinalResult.txt
perl line_select.pl $customized_data_path/target_customized.smi PDBID >target.smi
fi

mkdir tmpdir
mv * tmpdir/
cp tmpdir/FinalResult.txt .
cp tmpdir/target.smi .
date
fi
##EndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEndEnd
