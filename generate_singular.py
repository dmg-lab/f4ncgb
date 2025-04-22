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
    
    output.write('LIB "freegb.lib";\n')
    output.write('system("--ticks-per-sec", 100);\n')
    
    vars = '('
    for v in reversed(input.readline().rstrip().split(",")):
        vars += v + ','
    vars = vars[:-1]
    vars += ')'    
    char = input.readline().rstrip()
    
    ideal = ''
    for line in input:
        ideal += line.rstrip()
    ideal += ';'

    output.write('ring rr = ' + char + ',' + vars + ',Dp;\n')
    output.write('ring RR = freeAlgebra(rr,' + d + ');\n')
    output.write('ideal II = ' + ideal + '\n')
    output.write('timer = 0; ideal GG = twostd(II); timer;')
    output.write('quit;')
