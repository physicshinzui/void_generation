#!/bin/bash
#set -Ceu
set -eu

cat << EOS

Usage: 
    bash $0 [PATH_TO_DATA_DIR including PDBs] 

EOS

DATA_DIR=$1 #data_holo
suffix=pdb
#suffix=ent 

files=$(ls -1 ${DATA_DIR}/*.${suffix})
mkdir outputs

function make_void () {
    inp='templates/input_void_no'
    pdb=$1
    name=$(basename $pdb)
    cat $inp | sed -e "s!#{INPDB}!${pdb}!g" \
                   -e "s!#{OUTPDB}!outputs/${name%.*}_void.pdb!g" > outputs/${name%.*}_void.inp
    
    src/pdb_void_disp < outputs/${name%.*}_void.inp > outputs/${name%.*}.log
    mv void.pdb outputs/void_${name}
}

for file in $files; do
    echo "--$file--"
    time make_void $file
    echo "---------"
done

