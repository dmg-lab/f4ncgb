#!/usr/bin/env python3

import sys
import getopt

def read_input(argv):
 
    opts, args = getopt.getopt(argv[1:], "i:o:d:")
    
    input = None
    output = None
    deg = 0
        
    for opt, arg in opts:
        if opt == "-i":
            input = open(arg, 'r')
        elif opt == "-o":
            output = open(arg,'w')
        elif opt == "-d":
            deg = arg
        
    return input, output, deg
    
if __name__ == '__main__':  
    
    input, output, d = read_input(sys.argv)
    
    output.write('LIB "freegb.lib";\n')
    output.write('system("--ticks-per-sec", 100);\n')
    
    vars = '('
    for v in input.readline().rstrip().split(","):
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
    output.write('timer = 0; ideal G = twostd(I); timer;')
    output.write('quit;')