#!/usr/bin/env python3

import sys
import getopt

import argparse

def read_input(argv):
    parser = argparse.ArgumentParser(description='Process some inputs.')
    parser.add_argument('input', type=argparse.FileType('r'), help='Input file')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), required=True, help='Output file')
    parser.add_argument('-d', dest='deg', type=int, default=0, help='Degree value')

    args = parser.parse_args(argv[1:])
    
    return args.input, args.output, str(args.deg)
    
if __name__ == '__main__':  
    
    input, output, d = read_input(sys.argv)
    
    output.write('from operator_gb import *\n')
    
    vars = '['
    for v in input.readline().rstrip().split(","):
        vars += '"' + v + '",'
    vars = vars[:-1]
    vars += ']'    
    char = input.readline().rstrip()
    
    ideal = '['
    for line in input:
        ideal += 'F("' + line.rstrip().replace(',','') + '"), '
    ideal += ']'
    
    if char == '0':
        char = 'QQ'
    else:
        char = 'GF(' + char + ')'  
    
    output.write('F = FreeAlgebra(' + char + ',' + vars + ')\n')
    output.write('I = NCIdeal(' + ideal + ')\n')
    output.write('G = I.groebner_basis(10**5,maxdeg=' + d + ', trace_cofactors=False, verbose=1)\n')
    output.write('quit')
    
