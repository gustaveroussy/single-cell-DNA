#!/bin/bash

cd dir_bin/infSCITE

csv_file=/mnt/data/infscite_parameters.csv

ado=$(cat $csv_file | cut -d"," -f1 | tail -n1)
fdr=$(cat $csv_file | cut -d"," -f2 | tail -n1)
doublet_rate=$(cat $csv_file | cut -d"," -f3 | tail -n1)
n_sample=$(cat $csv_file | cut -d"," -f4 | tail -n1)
n_variants=$(cat $csv_file | cut -d"," -f5 | tail -n1)


./infSCITE -i /mnt/data/mutation_matrix.csv -n $n_variants -m $n_sample -r 1 -l 5000 -fd $fdr -ad $ado -s -e 0.2 -p 10000 -d 0.08 -o /mnt/output/$2 -names /mnt/data/gene_list.geneNames -seed 1