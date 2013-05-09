#!/bin/bash

#position|polymorphism|corr_coeff

corr_coeff=$1
ref_file=$2

num_coeffs=$(cat $ref_file | grep -v "Mutation corr_coeff" | wc -l)

grep -v "Mutation corr_coeff" $ref_file |
    cut -d' ' -f2 > tmp_$$.txt
echo $corr_coeff >> tmp_$$.txt

pos=$(sort -n < tmp_$$.txt | grep -n "$corr_coeff" | cut -d: -f 1)
pos=$(echo $pos | cut -d' ' -f1)
rm tmp_$$.txt
if [ $pos -gt $num_coeffs ]
then
	pos=$num_coeffs
fi

rank=$(echo "$pos/$num_coeffs")
p_value=$(echo "scale=4; $pos / $num_coeffs" | bc)

echo "$rank:$p_value"
