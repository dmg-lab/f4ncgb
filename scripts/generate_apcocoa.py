#!/usr/bin/env python3

import sys
import getopt
import re

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
    
    output.write('Alias NC := $apcocoa/ncpoly;\n')
    
    vars = '['
    for v in input.readline().rstrip().split(","):
        vars += v + ','
    vars = vars[:-1]
    vars += ']'   
    
    char = input.readline() 
    
    ideal = []
    for line in input:
        ideal.append(line.rstrip().replace(",",""))
    
    output.write('Use QQ' + vars + ';\n')
    output.write('NC.SetOrdering("LLEX");\n')
    i = 1
    for f in ideal:
        f = f.replace("+",",").replace("-", ",-")
        f = f.split(",")
        output.write('F' + str(i) + ' := [')
        first = True
        for m in f:
            if not first:
                output.write(',')
            first=False
            output.write('[')
            vv = ''
            for v in m.split("*"):
                vv += v + ','
            output.write(vv[:-1])
            output.write(']')
        output.write('];\n')
        i+=1
    
    I = ''
    for i in range(1,len(ideal)+1):
        I += 'F' + str(i) + ','
    output.write('G := [' + I[:-1] + '];\n')
    output.write('NC.TruncatedGB(G, ' + d + ', 31, 0);\n')
    output.write('Quit;\n')
