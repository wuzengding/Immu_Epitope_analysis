##########################################################################
# File Name:        TransMembrane.py
# Author:           wuzengding
# mail:             wuzengding@163.com
# Created Time:     Thu 12 Oct 2022 16:05:25 PM CST
##########################################################################

import time
import os, os.path
from argparse import ArgumentParser
import subprocess
import pandas as pd



if __name__ == '__main__':
    data_start = time.time()
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())) + " TransMembrane bed making .....")
    parser = ArgumentParser()
    parser.add_argument("-f","--TMRf",dest="transmembrane_file",
        help="file of bastp result",required=True)
    parser.add_argument("-o","--outdir",dest="outdir",
        help="out path",required=True)
    args = parser.parse_args()
    TMRf = args.transmembrane_file
    outdir = args.outdir
    
    with open(outdir+"/TransMembrane.bed","w") as tmrbed:
        tmrbed.write("""track name="TransMembrance" description="TransMembrance" visibility=2 itemRgb="On" useScore=1"""+"\n")
        with open(TMRf,"r") as tmrf:
            for line in tmrf:
                print(line)
                if not (line.startswith("#") or line.startswith("//")):
                    print(line)
                    qseqid,start,end,sseqid,strand,thickStart,thickEnd = \
                    line.split('\t')[0],line.split('\t')[2],line.split('\t')[3],line.split('\t')[1], \
                    ".",line.split('\t')[2],line.split('\t')[3]
                    if sseqid != "TMhelix":
                        score="100"
                        itemRgb="0,255,0"
                    else:
                        score="1000"
                        itemRgb="255,0,0"
                    tmrbed.write('\t'.join([qseqid,start,end,sseqid,score,strand,thickStart,thickEnd,itemRgb])+"\n")