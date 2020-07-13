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

#>  Author of the script: Goutam Mukherjee
#>  Purpose: Generate physicochemical parameter of ligand (drug-like) molecule
#>  How to run:
#>  bash RASPD_Ligand_parameter.sh <ligand id>
#>  Say,
#>  ligand id is: lig.pdb
#>  bash RASPD_Ligand_parameter.sh lig
#>  Please note that "lig.pdb" file must be present where the script, "RASPD_Ligand_parameter.sh" is executed.
#>  Output file is:
#>  "MOLECULARPROPERTY" contains an information about physicochemical properties of query ligand molecule (say, lig.pdb)

if [ -z "$raspd_root" ]
then
	echo -e '\033[1mPlease set path of the RASPD+ repository, miniconda and TRAPP to <path_to_RASPD+_repository>/config/init.sh file and source it\033[0m'
        echo -e '\033[1m   source <path_to_RASPD+_repository>/config/init.sh\033[0m'
        echo " "
        echo -e '\033[1m You should not run this script directly\033[0m'
        echo -e '\033[1m Please go to the ../copy/ folder and run the script lig_parameters_gen.sh\033[0m'
        echo -e '\033[1m lig_parameters_gen.sh call this script to calculate physicochemical properties of small molecules\033[0m'
	echo -e " "	
	echo -e '\033[1m However, you can run this script by the following command\033[0m'	
	echo -e "bash RASPD_Ligand_parameter.sh lig"
	echo -e "where lig is a 3D coordinates of small molecule in pdb format, i.e.; input file name is lig.pdb. Please note, .pdb extension is not required here while specifying the input"
        exit
fi

export path=$raspd_root/bin
export scripts_path=$raspd_root/scripts/chem

line=$1
$path/pdbtopdb.exe $line".pdb" $line"K.pdb"
$path/mass.exe $line"K.pdb" >>Molecular_Weight
$path/Connect2.0.exe $line"K.pdb" $line".con" 1 >R
sed -i 's/000   0  0  0/001   0  0  0/g' $line".con"
$path/windex_single.exe $line".con" >$line".wi"

sz=`wc -l $line".wi" | awk '{print($1)}'`
if [ $sz -eq 0 ]
then
egrep "ATOM|CONNECT|0  0  0" $line".con" >temp.con
sed -i 's/000   0  0  0/001   0  0  0/g' temp.con
$path/windex_single.exe temp.con >$line".wi"
fi

$path/pdbarrange.exe $line"K.pdb" arrangeline >$line"MD.pdb"
$path/pdbtopdb.exe $line"MD.pdb" $line"M.pdb"
$path/Connect2.0.exe $line"M.pdb" $line"New.pdb" 1 > RR
sed -i 's/000   0  0  0/001   0  0  0/g' $line"New.pdb"
$path/Connect2.0.exe $line"M.pdb" $line"New1.pdb" 0 > RR
sed -i 's/000   0  0  0/001   0  0  0/g' $line"New1.pdb"

$path/lp1.exe $line"M.pdb" $line".aa"
$path/lp2.exe $line"New.pdb" $line".aa" $line".bb"
$path/lp3.exe $line"New.pdb" $line".cc"
/usr/bin/perl $scripts_path/lp4.pl $line"New1.pdb" >bondorder
$path/lp5.exe bondorder $line".bb" $line".cc" $line".dd"
$path/lp6.exe $line".dd" $line".cc" $line".ee"
$path/lp7.exe $line".ee" LPMR1 >RR
$path/lp8.exe LPMR1 LOGP MR
$path/hbd_hba.exe $line".bb" hydrogen_bond_donor hydrogen_bond_acceptor
$path/lipi.exe hydrogen_bond_donor hydrogen_bond_acceptor LOGP MR  $line".wi" >>MOLECULARPROPERTY

rm R $line".con" $line".wi" $line"K.pdb" $line"MD.pdb" $line"M.pdb" $line"New.pdb" $line"New1.pdb" RR $line".aa" $line".bb" $line".cc" $line".dd" $line".ee" LPMR1
rm output.log bondorder arrangeline LOGP MR hydrogen_bond_donor hydrogen_bond_acceptor tmp*
