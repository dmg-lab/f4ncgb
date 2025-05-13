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
    
    output.write('LoadPackage( "GBNP" );;\n')
    output.write('SetInfoLevel(InfoGBNP,1);\n')
    
    vars = '['
    V = []
    for v in input.readline().rstrip().split(","):
        vars += '"' + v + '",'
        V.append(v)
    vars = vars[:-1]
    vars += ']'    
    char = input.readline().rstrip()
    
    ideal = []
    for line in input:
        ideal.append('GP2NP(' + line.rstrip().replace(',','') + '), ')
    
    output.write('F := FreeAssociativeAlgebraWithOne(Rationals,' + vars + ');;\n')
    for v in V:
        output.write(v + ' := F.' + v + ';;\n')
    output.write('I := [')
    for f in ideal:
        output.write(f);
    output.write('];;\n')
    order = str([1]*len(V))
    output.write('G := SGrobnerTrunc(I, ' + d  + ' , ' + order + ');\n')
    output.write('quit;')
    
