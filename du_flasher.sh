#!/bin/bash

ip=10.0
cmdr=./cmdr

du=$1
clbarray=(101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118)


for base in "${clbarray[@]}"
do

address=${ip}.${du}.${base}

#flashing the runtime
${cmdr} ${address} img.prog --pos 3 /home/km3net/applications/clbng/firmware/


done
