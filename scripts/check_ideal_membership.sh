#!/bin/bash

sage ./check_basis.py $1 -o input.sing -b $2 -d $3

singular input.sing

rm input.sing