#!/usr/bin/env python

import os
import sys

# replace +/- and fragment

def main(infile, outfile, temp):

	awk1 = '{OFS=" "} {print $1, $2, $3, $4, 0, $6, $7, $8, 1, 10, 10}'
	awk2 = '{if ($3 > $7){ print $1, $6, $7, $8, $9, $2, $3, $4, $5, $10, $11}else {print}}'

	cmd = f'cat {infile} | sed \'s/+/0/g\' | sed \'s/-/1/g\' | awk \'{awk1}\' | awk \'{awk2}\' | sort -k3,3d -k7,7d -T {temp} > {outfile}'

	print(cmd)
	os.system(cmd)

if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2],sys.argv[3])
