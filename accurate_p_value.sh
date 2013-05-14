#!/bin/bash

#	Copyright 2010-2013 Matt Halvorsen, Sam Broadaway, J.S. Martin, Chas Kissick
#	This file is part of SNPFold.
#
#	SNPFold is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	SNPFold is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with SNPFold.  If not, see <http://www.gnu.org/licenses/>.


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
