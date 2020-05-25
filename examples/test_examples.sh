#!/bin/bash

cd `dirname $0`

python ../consensus_creator.py input_example.txt > tmp.txt

diff tmp.txt expected_output.txt
rm tmp.txt
