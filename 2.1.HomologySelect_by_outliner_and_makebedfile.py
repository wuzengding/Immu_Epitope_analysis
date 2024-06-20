##########################################################################
# File Name:        HomologySeleclt_by_outliner_and_makebedfile.py
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
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())) + " Homology Selecting .....")
    parser = ArgumentParser()
    parser.add_argument("-f","--blastf",dest="blast_file",
        help="file of bastp result",required=True)
    parser.add_argument("-i","--identity",dest="identity",
        help="percent of identity, default=65",required=False,default=65)
    parser.add_argument("-e","--evalue",dest="evalue",
        help="evalue cutoff, default=0",required=False,default=0)
    parser.add_argument("-b","-bitscore",dest="bitscore",
        help="bit score cutoff,default=0",required=False,default=0)
    parser.add_argument("-o","--outdir",dest="outdir",
        help="out path",required=True)
        
    args = parser.parse_args()
    blast_file =  args.blast_file
    identity = int(args.identity)
    evalue = int(args.evalue)
    bitscore = int(args.bitscore)
    outdir = args.outdir
    
    
    def prepend_line(file_name, line):
        """ Insert given string as a new line at the beginning of a file """
        # define name of temporary dummy file
        dummy_file = file_name + '.bak'
        # open original file in read mode and dummy file in write mode
        with open(file_name, 'r') as read_obj, open(dummy_file, 'w') as write_obj:
            # Write given line to the dummy file
            write_obj.write(line + '\n')
            # Read lines from original file one by one and append them to the dummy file
            for line in read_obj:
                write_obj.write(line)
        # remove original file
        os.remove(file_name)
        # Rename dummy file as the original file
        os.rename(dummy_file, file_name)
    
    
    headname= ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"]
    blastdf = pd.read_csv(blast_file,sep="\t",names=headname)
    blastdf[["protein","pos"]] = blastdf.qseqid.str.split(":",expand=True) 
    proteinlist = blastdf.protein.unique().tolist()
    selfUniproList = []
    for proteinid in proteinlist:
        # calculate the self 
        qblastdf = blastdf[blastdf["protein"]==proteinid]
        homoRefseries = qblastdf.sseqid.value_counts()
        # calculate the query number (or bin numbers)
        querynum = len(qblastdf.qseqid.unique())
        # find those reference id which is the query ploymer 
        selfUniproID = homoRefseries[homoRefseries>=querynum*0.95].index.tolist()
        selfUniproList += selfUniproID

    # remove those aligned lines which reference is it's self
    noSelfBlastdf = blastdf[~blastdf.sseqid.isin(selfUniproList)]
    # remove those aligned lines which pidnet bigger than cutoff
    noSelfHomodf = noSelfBlastdf[noSelfBlastdf.pident > identity]

    noSelfHomodf[["start","end"]] = noSelfHomodf.pos.str.split("-",expand=True)
    noSelfHomodf.loc[:,"score"] = noSelfHomodf.pident*10
    noSelfHomodf.loc[:,"strand"] = "."
    noSelfHomodf.loc[:,"thickStart"] = noSelfHomodf.loc[:,"start"]
    noSelfHomodf.loc[:,"thickEnd"] = noSelfHomodf.loc[:,"end"]
    noSelfHomodf.loc[:,"itemRgb"] = "255,0,0"
    noSelfHomodf.loc[:,"score"] = noSelfHomodf["score"].apply(lambda x:int(x))
    noSelfHomodf[["qseqid","start","end","sseqid","score","strand","thickStart","thickEnd","itemRgb"]].to_csv(outdir+"/homo_peptide.bed",sep="\t",index=False,header=False)
    bedheader = '''track name="Homologous" description="Homologous with others proteins" visibility=2 itemRgb="On" useScore=1'''
    prepend_line(outdir+"/homo_peptide.bed", bedheader)
    