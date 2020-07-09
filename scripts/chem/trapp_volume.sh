#>  Copyright (c) 2020
#>  Heidelberg University and Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org)
#>  Schloss-Wolfsbrunnenweg 35
#>  69118 Heidelberg, Germany
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

##     Author: Jui-Hung Yuan, Heidelberg Institute for Theoretical Studies HITS gGmbH, and Goutam Mukherjee
##     This script calculates the volume of the ligand binding pocket of a protein

eval "$($conda_root/condabin/conda shell.bash hook)"

if [ -z "$raspd_root" ]
then
	echo -e '\033[1mPlease set path of the TRAPP to <path_to_RASPD+_repository>/config/init.sh file and source it\033[0m'
        echo -e '\033[1m   source <path_to_RASPD+_repository>/config/init.sh\033[0m'
        echo " "
        echo -e '\033[1m You should not run this script directly\033[0m'
        echo -e '\033[1m This programme call TRAPP to calculate pocket volume of a target protein\033[0m'
        exit
fi

# TODO: Make the same assertion for path to TRAPP installation

line=$1
cwd=`pwd`
rm -rf $line
mkdir $line
egrep -v "DRG" $line".pdb" >$line/$line"_amber.pdb"
egrep "DRG" $line".pdb" >$line/$line"_sdf.pdb"

echo "$cwd/$line/$line" >info.txt
conda activate trappenv
echo "$TRAPP/scripts/smallset_genX.py"
python $TRAPP/scripts/smallset_genX.py --info info.txt --out goutam
conda deactivate
awk '{print($2)}' goutam_volume.txt >volm
