#!/bin/bash
# transfer gcmt Mrtp to Mxx,Mxy,Mxz,Myy,Myz,Mzz
gawk -v m1="$1" -v m2="$2" -v m3="$3" -v m4="$4" -v m5="$5" -v m6="$6" 'BEGIN{print m2,-m6,m4,m3,-m5,m1}'
