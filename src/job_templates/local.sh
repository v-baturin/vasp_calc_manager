#!/bin/bash
source ~/.bashrc
nohup $GAUSS_EXEDIR/g16 < @INPUT > @OUTFILE 2>&1 &
echo $!

