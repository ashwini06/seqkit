#!/usr/bin/python

import sys

fl = open("/proj/b2012025/RAW_DATA/ChIP_Data/Ascl1_US/Mark_Mash1_s4/motif/denovo/Ascl1_Foxa2boundpositions.txt",'r')
op_fl = open("/proj/b2012025/RAW_DATA/ChIP_Data/Ascl1_US/Mark_Mash1_s4/motif/denovo/Ascl1_Foxa2boundpositions_summits.txt",'w')
for ln in iter(fl):
    l = ln.strip()
    pos = l.split(":")
    summit = (int(pos[1])+int(pos[2]))/2  
    op_fl.write('{}\n'.format(summit))

    

