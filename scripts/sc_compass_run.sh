#!/bin/bash

cd dir_bin/COMPASS/

./COMPASS -i $1 -o $2 --nchains $3 --chainlength $4 --CNV $5
dot -Tpng -o $6 $7