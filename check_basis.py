#!/usr/bin/env python3

import sys
import getopt

import argparse

from sage.all import *

def read_input(argv):
    parser = argparse.ArgumentParser(description='Process some inputs.')
    parser.add_argument('input', type=argparse.FileType('r'), help='Input file')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), required=True, help='Output file')
    parser.add_argument('-b', '--basis', type=argparse.FileType('r'), required=True, help='Basis file')
    parser.add_argument('-d', dest='deg', type=int, default=0, help='Degree value')

    args = parser.parse_args(argv[1:])
    
    return args.input, args.output, args.basis, str(args.deg)
    
if __name__ == '__main__':  
    
    input, output, basis, d = read_input(sys.argv)
    
    output.write('LIB "freegb.lib";\n')
    
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
    
    output.write('ring r = ' + char + ',' + vars + ',Dp;\n')
    output.write('ring R = freeAlgebra(r,' + d + ');\n')
    output.write('ideal I = ' + ideal + '\n')
    output.write('ideal G = twostd(I);')
    output.write('write(":w basis.sing", G);\n')
    output.write('quit;')
    
    GB_sing = None
    GB_us = None
    
    try:
        with open('basis.sing', 'r') as f:
            GB_sing = f.readline().split(",")
    except:
        sys.exit(1)
        
    K = QQ
    char = int(char)
    if char != 0:
        K = GF(char)
    vars = list(reversed(vars[1:-1].split(",")))
    R = FreeAlgebra(QQ,vars)
            
    GB_sing = [R(g.rstrip().strip(",")) for g in GB_sing]
    GB_us = [R(g.rstrip().strip(",")) for g in basis]
            
    GB_sing = [g.leading_item()[0].to_word() for g in GB_sing]
    GB_us = [g.leading_item()[0].to_word() for g in GB_us]
    
    for g in GB_sing:
        divisible = False
        for f in GB_us:
            if str(f) in str(g):
                divisible = True
                break
        if not divisible:
            print("VERIFICATION FAILED" , g, "not in ideal")
            sys.exit(1)
    
    print("VERIFICATION SUCCEEDED")
            
        


    
            