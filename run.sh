#!/bin/bash

#./tf-cluster -1 1.5 -2 1.2 -3 0.8 -c "spearman" \
#  -e correlation-matrix/example/rma_salt.txt \
#  -t correlation-matrix/example/genelist.txt -k 100

valgrind -v --leak-check=full ./tf-cluster -1 1.5 -2 1.2 -3 0.8 -c "spearman" -e At_stem_rma_july2012_AGI129_final.txt -t AT_TF_list.csv -k 100
