#!/bin/bash

gfortran -O3 tb1v11.f90 -llapack -o tb1v11.exe
sleep 2 # wait 2 seconds
nohup time ./tb1v11.exe & 
